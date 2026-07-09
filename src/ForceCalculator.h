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
    double          dist_um         = 0.0;

    // Calculated values
    double          live_L_um       = std::numeric_limits<double>::quiet_NaN();
    double          var_x_um2       = std::numeric_limits<double>::quiet_NaN();
    double          fluct_x_um      = std::numeric_limits<double>::quiet_NaN();
    double          fluct2_x_um2    = std::numeric_limits<double>::quiet_NaN();
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

        FrameData fd{frame_, timestamp_ms, x1_px, y1_px, x2_px, y2_px, dist_um};

        // Update rolling window
        window_.push_back(fd);

        const size_t TARGET_FRAMES = VARIANCE_WINDOW + 1;
        size_t center_idx = VARIANCE_WINDOW / 2;
        // First 300 frames NaN force
        if (frame_ <= center_idx) {
            appendRecord(computeRecord(fd));
        }

        if (window_.size() > TARGET_FRAMES) {
            window_.pop_front();
        }

        if (window_.size() == TARGET_FRAMES) {
            ForceRecord rec = computeRecord(window_[center_idx]);
            out_force_pN = rec.force_pN;
            out_live_L_um = rec.live_L_um;
            appendRecord(rec);
        }

        if (CSV_FLUSH_INTERVAL > 0 && frame_ % CSV_FLUSH_INTERVAL == 0) {
            flushCSV();
        }
    }

    void finalize() {
        // Tail flush
        while (!window_.empty()) {
            const auto& fd = window_.front();
            if (fd.frame > last_logged_frame_) {
                appendRecord(computeRecord(fd));
            }
            window_.pop_front();
        }
        flushCSV(true);
    }

    void setMicronsPerPixel(double mpp) { dpp_local_ = mpp; }

    void reset() {
        frame_ = 0;
        last_logged_frame_ = 0;
        window_.clear();
        pending_records_.clear();
        file_header_written_ = false;
        if (ofs_.is_open()) ofs_.close();
    }

private:
    ForceRecord computeRecord(const FrameData& target) const {
        ForceRecord rec;
        rec.frame = target.frame;
        rec.timestamp_ms = target.timestamp_ms;
        rec.x1_px = target.x1_px; rec.y1_px = target.y1_px;
        rec.x2_px = target.x2_px; rec.y2_px = target.y2_px;
        rec.dist_um = target.dist_um;

        int valid_n = 0;
        double sum_x2_px = 0.0, sum_dist_um = 0.0;

        for (const auto& w : window_) {
            if (!std::isnan(w.x2_px) && !std::isnan(w.dist_um)) {
                sum_x2_px += w.x2_px;
                sum_dist_um += w.dist_um;
                valid_n++;
            }
        }

        if (valid_n >= MIN_VARIANCE_SAMPLES) {
            double avg_x2_px = sum_x2_px / valid_n;
            double avg_dist_um = sum_dist_um / valid_n;

            double var_sum_m2 = 0.0;
            double mpp_m = dpp_local_ * 1e-6;

            for (const auto& w : window_) {
                if (!std::isnan(w.x2_px) && !std::isnan(w.dist_um)) {
                    double dx_m = (w.x2_px - avg_x2_px) * mpp_m;
                    var_sum_m2 += (dx_m * dx_m);
                }
            }

            // We use n population variance to match dynforce
            double variance_m2 = var_sum_m2 / valid_n;

            rec.fluct_x_um = (target.x2_px - avg_x2_px) * dpp_local_;
            rec.fluct2_x_um2 = rec.fluct_x_um * rec.fluct_x_um;
            rec.live_L_um = avg_dist_um;
            rec.var_x_um2 = variance_m2 * 1e12;

            if (variance_m2 > 1e-24) {
                rec.force_pN = 1e12 * (kT_ * (avg_dist_um * 1e-6)) / variance_m2;
            }
        }
        return rec;
    }

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
            ofs_ << "frame,timestamp_ms,x1_px,y1_px,x2_px,y2_px,dist_um,live_L_um,var_x_um2,fluct_x_um,fluct2_x_um2,force_pN\n";
            file_header_written_ = true;
        }

        // Precision scale
        ofs_ << std::setprecision(9);
        for (const auto& r : pending_records_) {
            // Standard coordinates
            ofs_ << r.frame << ',' << r.timestamp_ms << ','
                 << r.x1_px << ',' << r.y1_px << ','
                 << r.x2_px << ',' << r.y2_px << ','
                 << r.dist_um << ',';

            // 0 (interchangeable) if NaN, else write
            if (std::isnan(r.live_L_um)) ofs_ << "0,"; else ofs_ << r.live_L_um << ',';
            if (std::isnan(r.var_x_um2)) ofs_ << "0,"; else ofs_ << r.var_x_um2 << ',';
            if (std::isnan(r.fluct_x_um)) ofs_ << "0,"; else ofs_ << r.fluct_x_um << ',';
            if (std::isnan(r.fluct2_x_um2)) ofs_ << "0,"; else ofs_ << r.fluct2_x_um2 << ',';
            if (std::isnan(r.force_pN)) ofs_ << "0\n"; else ofs_ << r.force_pN << '\n';
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