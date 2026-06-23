//
// Created by Taro Porschke on 6/22/26.
//

#ifndef FORCECALCULATOR_H
#define FORCECALCULATOR_H

// Uses Equipartition Theorem for real-time force estimation

// F = ( Thermal Energy * Average Extension ) / Variance
// where Thermal Energy is the Temperature * Boltzmann Constant
// where Variance is the transverse average squared fluctuation (degrees of freedom)

#pragma once
#include <cmath>
#include <deque>
#include <numeric>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "BeadTrackerCore.h"

// ----- Constants
// 297 K is constant temperature
// k_B * T at 297 K in pN * um
// Newton to piconewton = 1e12
// Meter to micrometer = 1e6
// The constant's name is long, so it is abbreviated somewhat
// K_Boltzmann pN * um/K
constexpr double kB_pN_um_perK = 1.380649e-23 * 1e12 * 1e6; // = 1.380649e-5 pN * um/K

// ----- Tunable Constants
constexpr int VARIANCE_WINDOW = 600;
constexpr int MIN_VARIANCE_SAMPLES = 120;
constexpr double DEFAULT_TEMPERATURE_K = 297.0;

// ----- CSV Output
// Written on destroyTracker call
// One row per frame
// Columns: frame, timestamp_ms, x1_px, y1_px, x2_px, y2_px,
// dist_um, delta_x_um, mean_L_um, var_x_um2, force_pN, n_variance_samples
// Baseline frames will have NaN force
constexpr int CSV_FLUSH_INTERVAL = 500;

struct FrameData {
    double x1, y1, x2, y2, dist;
};

struct ForceRecord {
    uint64_t        frame           = 0;
    double          timestamp_ms    = 0.0;
    float           x1_px           = 0.f;
    float           y1_px           = 0.f;
    float           x2_px           = 0.f;
    float           y2_px           = 0.f;
    double          dist_um         = 0.0;
    double          delta_x_um      = 0.0;
    double          live_L_um       = 0.0;
    double          var_x_um2       = 0.0;
    double          force_pN        = std::numeric_limits<double>::quiet_NaN();
};

// ----- Program
class ForceCalculator {
    public:
    void init(const std::string& csv_path, double temperature_K = DEFAULT_TEMPERATURE_K) {
        csv_path_ = csv_path;
        temperature_K_ = temperature_K;
        kT_ = kB_pN_um_perK * temperature_K_;
        reset();
    }

    void update(float x1_px, float y1_px, float x2_px, float y2_px, double dist_um,
                double timestamp_ms, double& out_force_pN)
    {
        out_force_pN = std::numeric_limits<double>::quiet_NaN();
        frame_++;

        ForceRecord rec;
        rec.frame = frame_;
        rec.timestamp_ms = timestamp_ms;
        rec.x1_px = x1_px; rec.y1_px = y1_px;
        rec.x2_px = x2_px; rec.y2_px = y2_px;
        rec.dist_um = dist_um;

        // Convert pixels to um
        FrameData fd{};
        fd.x1 = x1_px * MICRONS_PER_PIXEL_local_;
        fd.y1 = y1_px * MICRONS_PER_PIXEL_local_;
        fd.x2 = x2_px * MICRONS_PER_PIXEL_local_;
        fd.y2 = y2_px * MICRONS_PER_PIXEL_local_;
        fd.dist = dist_um;

        // Update rolling window
        window_.push_back(fd);
        if (window_.size() > VARIANCE_WINDOW) {
            window_.pop_front();
        }

        int n = static_cast<int>(window_.size());

        // Only calculate force if buffer full enough
        if (n >= MIN_VARIANCE_SAMPLES) {

            // Calculate local moving averages
            double m_x1 = 0, m_y1 = 0, m_x2 = 0, m_y2 = 0, m_L = 0;
            for (const auto& w : window_) {
                m_x1 += w.x1; m_y1 += w.y1;
                m_x2 += w.x2; m_y2 += w.y2;
                m_L += w.dist;
            }
            m_x1 /= n; m_y1 /= n; m_x2 /= n; m_y2 /= n; m_L /= n;

            // Local axis, midpoint
            double dx = m_x2 - m_x1;
            double dy = m_y2 - m_y1;
            double len = std::hypot(dx, dy);
            if (len < 1e-9) len = 1.0;

            // Perpendicular axis
            double axis_perp_x = -dy/len;
            double axis_perp_y = dx/len;

            double mid_x = (m_x1 + m_x2) * 0.5;
            double mid_y = (m_y1 + m_y2) * 0.5;

            // Calculate variance
            double var_sum = 0.0;
            double latest_delta_x = 0.0;

            for (const auto& w : window_) {
                double curr_mid_x = (w.x1 + w.x2) * 0.5;
                double curr_mid_y = (w.y1 + w.y2) * 0.5;

                // Displacement from moving mean
                double dmx = curr_mid_x - mid_x;
                double dmy = curr_mid_y - mid_y;

                // Project onto transverse axis
                double proj = dmx * axis_perp_x + dmy * axis_perp_y;
                var_sum += proj * proj;
                latest_delta_x = proj;
            }

            double variance = var_sum / (n - 1);

            // Calculate final force
            if (variance > 1e-9) {
                double force = (kT_ * m_L) / variance;
                rec.force_pN = force;
                out_force_pN = force;
            }

            rec.delta_x_um = latest_delta_x;
            rec.live_L_um = m_L;
            rec.var_x_um2 = variance;

        }

        appendRecord(rec);

        if constexpr (CSV_FLUSH_INTERVAL > 0) {
            if (frame_ % CSV_FLUSH_INTERVAL == 0) {
                flushCSV();
            }
        }
    }

    void finalize() { flushCSV(true); }
    void setMicronsPerPixel(double mpp) { MICRONS_PER_PIXEL_local_ = mpp; }

    void reset() {
        frame_ = 0;
        window_.clear();
        pending_records_.clear();
        file_header_written_ = false;
        if (ofs_.is_open()) ofs_.close();
    }

private:
    void appendRecord(const ForceRecord& r) { pending_records_.push_back(r);}

    void flushCSV(bool close = false) {
        if (csv_path_.empty() || pending_records_.empty()) {
            if (close && ofs_.is_open()) ofs_.close();
            return;
        }
        if (!ofs_.is_open()) ofs_.open(csv_path_, std::ios::app);
        if (!ofs_.is_open()) return;

        if (!file_header_written_) {
            ofs_ << "frame,timestamp_ms,dist_um,live_L_um,var_x_um2,forcepN\n";
            file_header_written_ = true;
        }

        ofs_ << std::fixed << std::setprecision(6);
        for (const auto& r : pending_records_) {
            ofs_ << r.frame << ',' << r.timestamp_ms << ',' << r.dist_um << ','
                << r.live_L_um << ',' << r.var_x_um2 << ',';
            if (std::isnan(r.force_pN)) ofs_ << '\n';
            else ofs_ << r.force_pN << '\n';
        }

        pending_records_.clear();
        if (close) ofs_.close();
    }

    double temperature_K_ = DEFAULT_TEMPERATURE_K;
    double kT_ = kB_pN_um_perK * DEFAULT_TEMPERATURE_K;
    double MICRONS_PER_PIXEL_local_ = 0.0623;
    uint64_t frame_ = 0;
    std::deque<FrameData> window_;
    std::string csv_path_;
    std::ofstream ofs_;
    bool file_header_written_ = false;
    std::vector<ForceRecord> pending_records_;
};




#endif //FORCECALCULATOR_H