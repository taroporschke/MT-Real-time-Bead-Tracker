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
#include <limits>

#include "BeadTrackerCore.h"

// ----- Constants
// 297 K is constant temperature
// k_B * T at 297 K in pN * um
// Newton to piconewton = 1e12
// Meter to micrometer = 1e6
// The constant's name is long, so it is abbreviated somewhat
// K_Boltzmann pN * um/K
// Matched to DynForce.cpp
constexpr double kB_J_perK = 1.3806503e-23;

// ----- Tunable Constants
constexpr int VARIANCE_WINDOW = 600;
constexpr int MIN_VARIANCE_SAMPLES = 120;
constexpr double DEFAULT_TEMPERATURE_K = 297.0;

// ----- CSV Output
// Written on destroyTracker call
// One row per frame
// Columns: frame, x1_px, y1_px, x2_px, y2_px,
// dist_um, delta_x_um, mean_L_um, var_x_um2, force_pN
// Baseline frames will have NaN force
constexpr int CSV_FLUSH_INTERVAL = 500;

struct FrameData {
    uint64_t frame;
    double timestamp_ms;
    float x1_px, y1_px, x2_px, y2_px;
    double dist_um;
};

struct ForceRecord {
    uint64_t        frame           = 0;
    double          timestamp_ms    = 0.0;
    float           x1_px           = 0.f;
    float           y1_px           = 0.f;
    float           x2_px           = 0.f;
    float           y2_px           = 0.f;
    double          dist_um         = 0.0;  // instantaneous distance
    double          delta_x_um      = 0.0;
    double          live_L_um       = 0.0;  // window avg distance
    double          var_x_um2       = 0.0;
    double          force_pN        = std::numeric_limits<double>::quiet_NaN();
};

// ----- Program
class ForceCalculator {
    public:
    void init(const std::string& csv_path, double temperature_K = DEFAULT_TEMPERATURE_K) {
        csv_path_ = csv_path;
        temperature_K_ = temperature_K;
        kT_ = kB_J_perK * temperature_K_;
        reset();
    }

    void update(float x1_px, float y1_px, float x2_px, float y2_px, double dist_um,
        double timestamp_ms, double& out_force_pN, double& out_live_L_um)
    {
        out_force_pN = std::numeric_limits<double>::quiet_NaN();
        out_live_L_um = std::numeric_limits<double>::quiet_NaN();
        frame_++;

        FrameData fd{};
        fd.frame = frame_;
        fd.timestamp_ms = timestamp_ms;
        fd.x1_px = x1_px; fd.y1_px = y1_px;
        fd.x2_px = x2_px; fd.y2_px = y2_px;
        fd.dist_um = dist_um;

        // Update rolling window
        window_.push_back(fd);

        // First 300 frames NaN force
        int center_idx = VARIANCE_WINDOW / 2;
        if (frame_ <= center_idx) {
            ForceRecord rec;
            rec.frame = fd.frame;
            rec.x1_px = fd.x1_px;
            rec.y1_px = fd.y1_px;
            rec.x2_px = fd.x2_px;
            rec.y2_px = fd.y2_px;
            rec.dist_um = fd.dist_um;
            rec.live_L_um = 0.0;
            rec.var_x_um2 = 0.0;
            rec.force_pN = 0.0;
            appendRecord(rec);
        }
        // -------------------------

        const int TARGET_FRAMES = VARIANCE_WINDOW + 1;
        if (window_.size() > TARGET_FRAMES) {
            window_.pop_front();
        }

        int n = static_cast<int>(window_.size());

        if (n == TARGET_FRAMES) {

            const auto& center_frame = window_[center_idx];

            ForceRecord rec;
            rec.frame = center_frame.frame;
            rec.timestamp_ms = center_frame.timestamp_ms;
            rec.x1_px = center_frame.x1_px; rec.y1_px = center_frame.y1_px;
            rec.x2_px = center_frame.x2_px; rec.y2_px = center_frame.y2_px;


            // Calculate moving averages
            // Not averaged to match DynForce
            double sum_x2_px = 0.0;
            double sum_dist_um = 0.0;

            for (const auto& w : window_) {
                sum_x2_px += w.x2_px;
                sum_dist_um += w.dist_um;
            }

            double avg_x2_px = sum_x2_px / n;
            double avg_dist_um = sum_dist_um / n;

            // Calculate variance
            double var_sum_m2 = 0.0;
            double mpp_m = dpp_local_ * 1e-6;

            for (const auto& w : window_) {
                double dx_px = w.x2_px - avg_x2_px;
                double dx_m = dx_px * mpp_m;
                var_sum_m2 += (dx_m * dx_m);
            }

            double variance_m2 = var_sum_m2 / n;

            // Calculate force
            double dist_m = center_frame.dist_um * 1e-6;
            double force_pN = std::numeric_limits<double>::quiet_NaN();

            if (variance_m2 > 1e-24) {
                // F = 1e12 * (KBT * dist) / variance
                force_pN = 1e12 * (kT_ * dist_m) / variance_m2;
            }

            out_force_pN = force_pN;
            out_live_L_um = avg_dist_um;

            // Record for CSV
            rec.dist_um = center_frame.dist_um;
            rec.live_L_um = avg_dist_um;
            rec.var_x_um2 = variance_m2 * 1e12;
            rec.force_pN = force_pN;

            // Frame's latest delta
            rec.delta_x_um = (center_frame.x2_px - avg_x2_px) * dpp_local_;

            appendRecord(rec);
        }

        if constexpr (CSV_FLUSH_INTERVAL > 0) {
            if (frame_ % CSV_FLUSH_INTERVAL == 0) {
                flushCSV();
            }
        }
    }

    void finalize() {
        // Catch remaining frames
        for (const auto& fd : window_) {
            if (fd.frame > last_logged_frame_) {
                ForceRecord rec;
                rec.frame = fd.frame;
                rec.x1_px = fd.x1_px;
                rec.y1_px = fd.y1_px;
                rec.x2_px = fd.x2_px;
                rec.y2_px = fd.y2_px;
                rec.dist_um = fd.dist_um;
                rec.live_L_um = 0.0;
                rec.var_x_um2 = 0.0;
                rec.force_pN = 0.0;
                appendRecord(rec);
            }
        }
        // ----------------------
        flushCSV(true);
    }

    void setMicronsPerPixel(double mpp) { dpp_local_ = mpp; }

    void reset() {
        frame_ = 0;
        last_logged_frame_ = 0; // Reset tracker memory
        window_.clear();
        pending_records_.clear();
        file_header_written_ = false;
        if (ofs_.is_open()) ofs_.close();
    }

private:
    void appendRecord(const ForceRecord& r) {
        pending_records_.push_back(r);
        // Keep track of highest frame queued
        if (r.frame > last_logged_frame_) {
            last_logged_frame_ = r.frame;
        }
    }

    void flushCSV(bool close = false) {
        if (csv_path_.empty() || pending_records_.empty()) {
            if (close && ofs_.is_open()) ofs_.close();
            return;
        }
        if (!ofs_.is_open()) ofs_.open(csv_path_, std::ios::app);
        if (!ofs_.is_open()) return;

        if (!file_header_written_) {
            ofs_ << "frame,x1_px,y1_px,x2_px,y2_px,dist_um,live_L_um,var_x_um2,forcepN\n";
            file_header_written_ = true;
        }

        ofs_ << std::fixed << std::setprecision(6);
        for (const auto& r : pending_records_) {
            ofs_    << r.frame << ','
                    << r.x1_px << ',' << r.y1_px << ','
                    << r.x2_px << ',' << r.y2_px << ','
                    << r.dist_um << ','
                    << r.live_L_um << ',' << r.var_x_um2 << ',';
            if (std::isnan(r.force_pN)) ofs_ << '\n';
            else ofs_ << r.force_pN << '\n';
        }

        pending_records_.clear();
        if (close) ofs_.close();
    }

    double temperature_K_ = DEFAULT_TEMPERATURE_K;
    double kT_ = kB_J_perK * DEFAULT_TEMPERATURE_K;
    double dpp_local_ = 0.0609;
    uint64_t frame_ = 0;
    uint64_t last_logged_frame_ = 0;
    std::deque<FrameData> window_;
    std::string csv_path_;
    std::ofstream ofs_;
    bool file_header_written_ = false;
    std::vector<ForceRecord> pending_records_;
};

#endif //FORCECALCULATOR_H