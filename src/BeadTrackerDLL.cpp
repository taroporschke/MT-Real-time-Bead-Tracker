//
// Created by Taro Porschke on 6/5/26.
//
#include "BeadTrackerDLL.h"
#include "BeadTrackerCore.h"
#include "ForceCalculator.h"
#include <opencv2/opencv.hpp>
#include <chrono>

// Internal state object - lives on heap between LabVIEW calls
// destroyTracker must be called to free
struct TrackerState {
    Bead bead1;
    Bead bead2;
    int width = 0;
    int height = 0;
    double lastFrameMs = 0.0;

    // --- New Engine Additions ---
    ForceCalculator force_eng;
    std::chrono::steady_clock::time_point startTime; // Tracks absolute experiment time
};

// ----- Create Tracker
extern "C" DLL_API long long createTracker(
    float x1, float y1,
    float x2, float y2,
    int width, int height,
    const char* csv_path,     // <-- Added CSV path parameter
    double temperature_K)     // <-- Added Temperature parameter
{
    try {
        // destroyTracker() calls delete on this pointer
        TrackerState* state = new TrackerState();

        state->width = width;
        state->height = height;

        // Start the experiment timer
        state->startTime = std::chrono::steady_clock::now();

        // Initialize bead 1
        state->bead1.center = cv::Point2f(x1, y1);
        state->bead1.lostFrames = 0;
        state->bead1.roi = clampROI(
            cv::Rect(
                (int)(x1 - BASE_ROI_TOP / 2.0),
                (int)(y1 - BASE_ROI_TOP / 2.0),
                BASE_ROI_TOP, BASE_ROI_TOP),
            cv::Size(width, height)
        );

        // Initialize bead 2
        state->bead2.center = cv::Point2f(x2, y2);
        state->bead2.lostFrames = 0;
        state->bead2.roi = clampROI(
            cv::Rect(
                (int)(x2 - BASE_ROI_BOTTOM / 2),
                (int)(y2 - BASE_ROI_BOTTOM / 2),
                BASE_ROI_BOTTOM, BASE_ROI_BOTTOM),
            cv::Size(width, height)
        );

        // Initialize ForceCalculator
        std::string path_str = (csv_path != nullptr) ? std::string(csv_path) : "";
        state->force_eng.init(path_str, temperature_K);

        // Ownership transferred to caller via opaque handle
        return reinterpret_cast<long long>(state);
    }
    catch (...) {
        return 0;
    }
}

// ----- TrackFrame
extern "C" DLL_API void trackFrame(
    long long      handle,
    unsigned char* pixels,
    int            width,
    int            height,
    float* out_x1,
    float* out_y1,
    float* out_x2,
    float* out_y2,
    float* out_dist_px,
    float* out_dist_um,
    double* out_force_pN)
{
    if (handle == 0 || pixels == nullptr) return;

    TrackerState* state = reinterpret_cast<TrackerState*>(handle);

    // Start timing for frame processing speed
    auto t1 = std::chrono::high_resolution_clock::now();

    // Calculate exact timestamp (in milliseconds) since createTracker was called
    double timestamp_ms = std::chrono::duration<double, std::milli>(t1 - state->startTime).count();

    cv::Mat gray(height, width, CV_8UC1, pixels);
    cv::Mat contourDisp = cv::Mat::zeros(gray.size(), CV_8UC3);

    // --- Track bead 1
    auto [c1, f1] = refineCenter(
        gray,
        state->bead1.roi,
        state->bead1.center,
        contourDisp,
        cv::Scalar(0, 0, 255)
    );

    // --- Track bead 2
    auto [c2, f2] = refineCenter(
        gray,
        state->bead2.roi,
        state->bead2.center,
        contourDisp,
        cv::Scalar(0, 255, 0)
    );

    // --- Update bead 1 state
    if (f1) {
        state->bead1.center = c1;
        state->bead1.lostFrames = 0;
        if (state->bead1.templ.empty()) {
            cv::Rect tR(static_cast<int>(c1.x - 15), static_cast<int>(c1.y - 15), 30, 30);
            tR &= cv::Rect(0, 0, gray.cols, gray.rows);
            state->bead1.templ = gray(tR).clone();
        }
        state->bead1.roi = clampROI(
            cv::Rect(
                (int)(c1.x - BASE_ROI_TOP / 2),
                (int)(c1.y - BASE_ROI_TOP / 2),
                BASE_ROI_TOP, BASE_ROI_TOP),
            gray.size()
        );
    }
    else {
        state->bead1.lostFrames++;
        if (state->bead1.lostFrames > 10)
            state->bead1.roi = cv::Rect(0, 0, width, height);
    }

    // --- Update bead 2 state
    if (f2) {
        state->bead2.center = c2;
        state->bead2.lostFrames = 0;
        if (state->bead2.templ.empty()) {
            cv::Rect tR(static_cast<int>(c2.x - 15), static_cast<int>(c2.y - 15), 30, 30);
            tR &= cv::Rect(0, 0, gray.cols, gray.rows);
            state->bead2.templ = gray(tR).clone();
        }
        state->bead2.roi = clampROI(
            cv::Rect(
                (int)(c2.x - BASE_ROI_BOTTOM / 2),
                (int)(c2.y - BASE_ROI_BOTTOM / 2),
                BASE_ROI_BOTTOM, BASE_ROI_BOTTOM),
            gray.size()
        );
    }
    else {
        state->bead2.lostFrames++;
        if (state->bead2.lostFrames > 10)
            state->bead2.roi = cv::Rect(0, 0, width, height);
    }

    // --- Compute distance
    double dist_px = cv::norm(state->bead1.center - state->bead2.center);
    double dist_um = dist_px * MICRONS_PER_PIXEL - BEAD_OFFSET_UM;

    // --- Feed physics engine & get Force
    double force_pN = 0.0;
    state->force_eng.update(
        state->bead1.center.x, state->bead1.center.y,
        state->bead2.center.x, state->bead2.center.y,
        dist_um, timestamp_ms, force_pN
    );

    // --- Stop timing
    auto t2 = std::chrono::high_resolution_clock::now();
    state->lastFrameMs = std::chrono::duration<double, std::milli>(t2 - t1).count();

    // --- Write results to LabVIEW output pointers (safely handle nullptrs)
    if (out_x1) *out_x1 = state->bead1.center.x;
    if (out_y1) *out_y1 = state->bead1.center.y;
    if (out_x2) *out_x2 = state->bead2.center.x;
    if (out_y2) *out_y2 = state->bead2.center.y;
    if (out_dist_px) *out_dist_px = (float)dist_px;
    if (out_dist_um) *out_dist_um = (float)dist_um;
    if (out_force_pN) *out_force_pN = force_pN;
}

// ----- destroyTracker
extern "C" DLL_API void destroyTracker(long long handle)
{
    if (handle == 0) return;
    TrackerState* state = reinterpret_cast<TrackerState*>(handle);

    // Flush remaining records to CSV before shutting down
    state->force_eng.finalize();

    delete state;
}

// ----- getLastFrameTimeMs
extern "C" DLL_API double getLastFrameTimeMs(long long handle)
{
    if (handle == 0) return -1.0;
    TrackerState* state = reinterpret_cast<TrackerState*>(handle);
    return state->lastFrameMs;
}