//
// Created by Taro Porschke on 6/4/26.
//
#include "BeadTrackerCore.h"
#include "BeadTrackerDLL.h"
#include <opencv2/opencv.hpp>
#include <iostream>
#include <chrono>

int main() {
    // Load test frame
    cv::Mat frame = cv::imread("C:/Users/rajendrab/Pictures/testframevlc.png");
    if (frame.empty()) {
        std::cerr << "Could not load image — check path\n";
        std::cin.get();
        return -1;
    }

    // Dimensions
    int width  = frame.cols;
    int height = frame.rows;
    std::cout << "Frame size: " << width << "x" << height << "\n";

    // IF NEEDED - convert to grayscale
    cv::Mat gray;
    if (frame.channels() == 3)
        cv::cvtColor(frame, gray, cv::COLOR_BGR2GRAY);
    else
        gray = frame;

    // Ensures continuity
    if (!gray.isContinuous())
        gray = gray.clone();

    // Hardcoded bead starting positions based on referenced test frame
    float x1 = 142.0f, y1 = 98.0f;
    float x2 = 149.0f, y2 = 481.0f;

    // Test 1 - Create Tracker
    std::cout << "\n--- TEST 1: createTracker ---\n";
    // Add all tests
    long long handle = createTracker(x1, y1, x2, y2, width, height, "sandbox_test_data.csv", 297.0);

    if (handle == 0) {
        std::cerr << "FAILED: createTracker returned 0\n";
        return -1;
    }
    std::cout << "PASSED: handle = " << handle << "\n";

    // Test 2 - Track Frame
    std::cout << "\n--- TEST 2: trackFrame (single frame) ---\n";

    float out_x1, out_y1, out_x2, out_y2, out_dist_px, out_dist_um;
    double force_pN = 0.0;

    trackFrame(
        handle,
        gray.data,         // raw pixel pointer
        width, height,
        &out_x1, &out_y1,
        &out_x2, &out_y2,
        &out_dist_px,
        &out_dist_um,
        &force_pN
    );

    std::cout << "Bead 1: (" << out_x1 << ", " << out_y1 << ")\n";
    std::cout << "Bead 2: (" << out_x2 << ", " << out_y2 << ")\n";
    std::cout << "Distance: " << out_dist_px << " px  |  "
              << out_dist_um << " um\n";
    std::cout << "Force: " << force_pN << " pN\n";

    // Sanity checks
    bool bead1_sane = (out_x1 > 0 && out_x1 < width &&
                       out_y1 > 0 && out_y1 < height);
    bool bead2_sane = (out_x2 > 0 && out_x2 < width &&
                       out_y2 > 0 && out_y2 < height);
    bool dist_sane  = (out_dist_px > 0 && out_dist_um > 0);

    std::cout << "Bead 1 position sane: " << (bead1_sane ? "YES" : "NO") << "\n";
    std::cout << "Bead 2 position sane: " << (bead2_sane ? "YES" : "NO") << "\n";
    std::cout << "Distance sane:        " << (dist_sane  ? "YES" : "NO") << "\n";

    // Test 3 - Last Frame Time
    std::cout << "\n--- TEST 3: getLastFrameTimeMs ---\n";
    double last_ms = getLastFrameTimeMs(handle);
    std::cout << "Last frame time: " << last_ms << " ms\n";
    std::cout << "Equivalent FPS: " << 1000.0 / last_ms << "\n";

    // Test 4 - 100 Frame Loop Timing
    std::cout << "\n--- TEST 4: 120-frame loop timing ---\n";

    auto t1 = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < 120; i++) {
        trackFrame(
            handle,
            gray.data,
            width, height,
            &out_x1, &out_y1,
            &out_x2, &out_y2,
            &out_dist_px,
            &out_dist_um,
            &force_pN
        );
    }

    auto t2 = std::chrono::high_resolution_clock::now();
    double total_ms = std::chrono::duration<double, std::milli>(t2 - t1).count();
    double per_frame_ms = total_ms / 120.0;
    double equiv_fps = 1000.0 / per_frame_ms;

    std::cout << "120 frames in " << total_ms << " ms\n";
    std::cout << "Per frame: " << per_frame_ms << " ms\n";
    std::cout << "Equivalent FPS ceiling: " << equiv_fps << "\n";
    std::cout << "8.3ms budget (120fps): "
              << (per_frame_ms < 8.3 ? "WITHIN BUDGET" : "OVER BUDGET") << "\n";

    // Test 5 - Destroy Tracker
    std::cout << "\n--- TEST 5: destroyTracker ---\n";
    destroyTracker(handle);
    std::cout << "PASSED: no crash\n";

    std::cout << "\n=== All tests complete ===\n";

    std::cout << "\n=== All tests complete ===\n";

    // Force console to wait for user
    std::cin.get();

    return 0;
}