//
// Created by Taro Porschke on 6/4/26.
//

#include "BeadTrackerCore.h"

// clampROI implementation
// Both corners are independently clamped, then w/h is derived
cv::Rect clampROI(const cv::Rect& r, const cv::Size& imsz) {
    // Clamp the top-left corner
    int x = std::max(0, r.x);
    int y = std::max(0, r.y);
    // Clamp the bottom-right corner
    int x2 = std::min(r.x + r.width, imsz.width);
    int y2 = std::min(r.y + r.height, imsz.height);
    // Derive width/height
    return { x, y, std::max(x2 - x, 0), std::max(y2 - y, 0) };
}

// refineCenter implementation
std::pair<cv::Point2f, bool> refineCenter(const cv::Mat& gray,
    cv::Rect roi,
    const cv::Point2f& prevCenter,
    cv::Mat& contourDisp,
    const cv::Scalar& color,
    int maxStepPx,
    double initialThreshFactor,
    const cv::Mat& beadTemplate)
{
    // --- Initialization ---
    cv::Rect imgRect(0, 0, gray.cols, gray.rows);
    cv::Point2f best = prevCenter; // Default to previous location
    bool ok = false;               // Success flag
    int expansion = 0;             // ROI expansion counter

    // Try multiple ROI expansions if the bead isn't found immediately
    while (!ok && expansion < 3) {
        int ex = std::max(0, roi.x - expansion * roi.width / 2);
        int ey = std::max(0, roi.y - expansion * roi.height / 2);
        int ex2 = std::min(roi.x + roi.width * (1 + expansion), gray.cols);
        int ey2 = std::min(roi.y + roi.height * (1 + expansion), gray.rows);

        cv::Rect expandedROI(ex, ey, ex2 - ex, ey2 - ey);
        expandedROI &= imgRect;

        // Skip if ROI became too small
        if (expandedROI.width <= 5 || expandedROI.height <= 5) break;

        // Preprocessing
        cv::Mat roiImg = gray(expandedROI);
        cv::Ptr<cv::CLAHE> clahe = cv::createCLAHE(3.0, cv::Size(8, 8)); // Contrast enhancement
        cv::Mat enhanced; clahe->apply(roiImg, enhanced);
        cv::Mat blurImg; cv::GaussianBlur(enhanced, blurImg, cv::Size(5, 5), 1.2); // Noise reduction
        cv::Mat norm8; cv::normalize(blurImg, norm8, 0, 255, cv::NORM_MINMAX, CV_8U); // Normalize to 8-bit

        // Adaptive threshold search
        double threshFactor = initialThreshFactor;
        const double minThreshFactor = 0.5;
        const double step = 0.1;

        while (threshFactor >= minThreshFactor && !ok) {
            int highThr = std::max(160, (int)std::round(threshFactor * 255));

            // Binary threshold and morphology cleanup
            cv::Mat bin;
            cv::threshold(norm8, bin, highThr, 255, cv::THRESH_BINARY);
            cv::morphologyEx(bin, bin, cv::MORPH_OPEN,
                cv::getStructuringElement(cv::MORPH_ELLIPSE, { 3,3 }));

            // Extract contours
            std::vector<std::vector<cv::Point>> cc;
            cv::findContours(bin, cc, cv::RETR_EXTERNAL, cv::CHAIN_APPROX_SIMPLE);

            // Lambda function to compute contour centroid
            auto centroidOf = [](const std::vector<cv::Point>& c)->cv::Point2f {
                cv::Moments m = cv::moments(c);
                return (m.m00 > 0) ?
                    cv::Point2f((float)(m.m10 / m.m00), (float)(m.m01 / m.m00)) :
                    cv::Point2f(-1, -1);
                };

            // Choose contour closest to previous center
            if (!cc.empty()) {
                double bestD = 1e12;
                for (auto& c : cc) {
                    double a = cv::contourArea(c);
                    if (a < 3 || a > 200) continue;   // Reject noise
                    cv::Point2f loc = centroidOf(c);
                    if (loc.x < 0) continue;
                    cv::Point2f cand(expandedROI.x + loc.x, expandedROI.y + loc.y);
                    double d = cv::norm(cand - prevCenter);
                    if (d < bestD) { bestD = d; best = cand; ok = true; }
                }
                // Discard if jump too large (to avoid false jumps)
                if (ok && cv::norm(best - prevCenter) > maxStepPx * (1 + expansion))
                    ok = false;
            }
            if (!ok) threshFactor -= step; // Lower threshold gradually
        }
        expansion++; // Increase ROI for next iteration if not found
    }

    // If contours fail, template matching
    if (!ok && !beadTemplate.empty()) {
        cv::Mat result;
        cv::matchTemplate(gray, beadTemplate, result, cv::TM_CCOEFF_NORMED);
        double minVal, maxVal;
        cv::Point minLoc, maxLoc;
        cv::minMaxLoc(result, &minVal, &maxVal, &minLoc, &maxLoc);
        if (maxVal > 0.5) {
            best = cv::Point2f(maxLoc.x + beadTemplate.cols / 2.0f,
                maxLoc.y + beadTemplate.rows / 2.0f);
            ok = true;
        }
    }

    return { best, ok };  // Return refined center and success flag
}