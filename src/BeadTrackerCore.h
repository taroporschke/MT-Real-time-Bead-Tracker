//
// Created by Taro Porschke on 6/4/26.
//


#ifndef BEADTRACKERCORE_H
#define BEADTRACKERCORE_H

#include <opencv2/opencv.hpp>

constexpr double MICRONS_PER_PIXEL = 0.0623;
constexpr double BEAD_OFFSET_UM = 2.8;
constexpr int BASE_ROI_TOP = 40;
constexpr int BASE_ROI_BOTTOM = 40;

struct Bead {
    cv::Point2f center;   // Centroid position (floating-point pixel coordinates)
    cv::Rect roi;         // Current ROI bounding box
    int lostFrames = 0;   // Counter for consecutive detection failures
    cv::Mat templ;        // Bead template (small image patch for fallback matching)
};

//clampROI declaration
cv::Rect clampROI(const cv::Rect& r, const cv::Size& imsz);

//refine center declaration
std::pair<cv::Point2f, bool> refineCenter(
    const cv::Mat& gray,
    cv::Rect roi,
    const cv::Point2f& prevCenter,
    cv::Mat& contourDisp,
    const cv::Scalar& color,
    int maxStepPx = 8,
    double initialThreshFactor = 0.85,
    const cv::Mat& beadTemplate = cv::Mat()
    );
#endif //BEADTRACKERCORE_H