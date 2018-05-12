/*#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/core/core.hpp"

using namespace std;
using namespace cv;
//#include "Functions.h"
int main()
{
Mat image;
//image = imread("lena.jpg",1);
image = imread("fp.tif",1);
if(image.empty())
{
cout << "Could not open or find the image" << std::endl ;
return -1;
}

/// Convert it to gray
//cvtColor( image, image, CV_RGB2GRAY );
//resize(image,image,Size(0,0),0.5,0.5,INTER_LINEAR);
namedWindow("Image", CV_WINDOW_AUTOSIZE );
imshow("Image", image);
/// Generate grad_x and grad_y
Mat grad_x, grad_y;
Mat abs_grad_x, abs_grad_y;
int scale = 1;
int delta = 0;
int ddepth = CV_16S;
Mat grad;


/// Gradient X
//Scharr( image, grad_x, ddepth, 1, 0, scale, delta, BORDER_DEFAULT );
Sobel( image, grad_x, ddepth, 1, 0, 3, scale, delta, BORDER_DEFAULT );
convertScaleAbs( grad_x, abs_grad_x );

/// Gradient Y
// Scharr( image, grad_y, ddepth, 0, 1, scale, delta, BORDER_DEFAULT );
Sobel( image, grad_y, ddepth, 0, 1, 3, scale, delta, BORDER_DEFAULT );
convertScaleAbs( grad_y, abs_grad_y );
/// Total Gradient (approximate)

addWeighted( abs_grad_x, 0.5, abs_grad_y, 0.5, 0, grad );




 // Mat orientation = Mat(abs_grad_x.rows, abs_grad_y.cols, CV_32F); //to store the gradients
  Mat img=Mat(abs_grad_x.rows, abs_grad_y.cols, CV_32F);//to draw out the map
  img = cv::Scalar(255,255,255);//all white


 // Calculate orientations of gradients --> in degrees
// Loop over all matrix values and calculate the accompanied orientation

 // for(int i = 0; i < abs_grad_x.rows; i++){
 //    for(int j = 0; j < abs_grad_x.cols; j++){
 //        // Retrieve a single value

 //        float valueX = abs_grad_x.at<float>(i,j);
 //        float valueY = abs_grad_x.at<float>(i,j);
 //        // Calculate the corresponding single direction, done by applying the arctangens function
 //        float result = fastAtan2(valueX,valueY);
 //        // Store in orientation matrix element
 //        orientation.at<float>(i,j) = result;

 //    }
 // }

  Mat orientation = Mat::zeros(abs_grad_x.rows, abs_grad_y.cols, CV_32F); //to store the gradients grad_x.convertTo(grad_x,CV_32F);
  grad_x.convertTo(grad_x,CV_32F);

grad_y.convertTo(grad_y,CV_32F);

phase(grad_x, grad_y, orientation, true);

cv::normalize(orientation, orientation, 0x00, 0xFF, cv::NORM_MINMAX, CV_8U);

namedWindow("Orientation", CV_WINDOW_AUTOSIZE );

imshow( "Orientation", orientation );

namedWindow("ImageSobel", CV_WINDOW_AUTOSIZE );
imshow( "ImageSobel", grad );

namedWindow("ImageSobelGx", CV_WINDOW_AUTOSIZE );
imshow( "ImageSobelGx", abs_grad_x );

namedWindow("ImageSobelGy", CV_WINDOW_AUTOSIZE );
imshow( "ImageSobelGy", abs_grad_y );

cout << orientation << endl;

waitKey(0);
return 0;
}
*/

#include <opencv2/opencv.hpp>
#include <iostream>
using namespace cv;
using namespace std;

#define THRESHOLD 10
#define block_len 16
#define block_width 16
#define pi 3.14
#define lowPassSize 5

int mask[256/block_width*256/block_len];

// function for splitting image into multiple blocks. rowDivisor and colDivisor specify the number of blocks in rows and cols respectively
int subdivide(cv::Mat &img, const int rowDivisor, const int colDivisor, std::vector<cv::Mat> &blocks) {
    /* Checking if the image was passed correctly */
    if(!img.data || img.empty())
        std::cerr << "Problem Loading Image" << std::endl;

    /* Cloning the image to another for visualization later, if you do not want to visualize the result just comment every line related to visualization */
    cv::Mat maskImg = img.clone();
    /* Checking if the clone image was cloned correctly */
    if(!maskImg.data || maskImg.empty())
        std::cerr << "Problem Loading Image" << std::endl;

    // check if divisors fit to image dimensions
    if(img.cols % colDivisor == 0 && img.rows % rowDivisor == 0) {
        // Iterator for mask
        int k =0;
        for(int y = 0; y < img.cols; y += img.cols / colDivisor) {
            for(int x = 0; x < img.rows; x += img.rows / rowDivisor) {
                blocks.push_back(img(cv::Rect(y, x, (img.cols / colDivisor), (img.rows / rowDivisor))));
                Scalar mean, stddev;
                meanStdDev(blocks[k], mean, stddev);

                // Draw red rectangle to show fingerprint area
                if (stddev.val[0] > THRESHOLD)
                    rectangle(maskImg, Point(y, x), Point(y + (maskImg.cols / colDivisor) - 1, x + (maskImg.rows / rowDivisor) - 1), CV_RGB(255, 0, 0), 1); // visualization
                k++;
            }
        }
        imshow("Segmented", maskImg); // visualization
        waitKey(0); // visualization
    } else if(img.cols % colDivisor != 0) {
        cerr << "Error: Please use another divisor for the column split." << endl;
        exit(1);
    } else if(img.rows % rowDivisor != 0) {
        cerr << "Error: Please use another divisor for the row split." << endl;
        exit(1);
    }
    return EXIT_SUCCESS;
}

int segment(Mat &block, float thres) {
    Scalar mean, stddev;
    meanStdDev(block, mean, stddev);
    if (stddev.val[0] < thres) {
        // block.setTo(Scalar(0,0,0));
        return 0;
    }
    return 1;
}

cv::Mat orientation(cv::Mat inputImage)
{
    cv::Mat orientationMat = cv::Mat::zeros(inputImage.size(), CV_8UC1);

    // compute gradients at each pixel
    cv::Mat grad_x, grad_y;
    cv::Sobel(inputImage, grad_x, CV_32F, 1, 0, 3, 1, 0, cv::BORDER_DEFAULT);
    cv::Sobel(inputImage, grad_y, CV_32F, 0, 1, 3, 1, 0, cv::BORDER_DEFAULT);

    cv::Mat Vx, Vy, theta, lowPassX, lowPassY;
    cv::Mat lowPassX2, lowPassY2;

    Vx = cv::Mat::zeros(inputImage.size(), inputImage.type());
    Vx.copyTo(Vy);
    Vx.copyTo(theta);
    Vx.copyTo(lowPassX);
    Vx.copyTo(lowPassY);
    Vx.copyTo(lowPassX2);
    Vx.copyTo(lowPassY2);

    cout << grad_x << endl;

    // estimate the local orientation of each block
    int blockSize = 16;

    for(int i = blockSize/2; i < inputImage.rows - blockSize/2; i+=blockSize)
    {    
        for(int j = blockSize / 2; j < inputImage.cols - blockSize/2; j+= blockSize)
        {
            float sum1 = 0.0;
            float sum2 = 0.0;

            for ( int u = i - blockSize/2; u < i + blockSize/2; u++)
            {
                for( int v = j - blockSize/2; v < j+blockSize/2; v++)
                {
                    sum1 += 2*grad_x.at<float>(u,v) * grad_y.at<float>(u,v);
                    sum2 += (grad_x.at<float>(u,v)*grad_x.at<float>(u,v)) - (grad_y.at<float>(u,v)*grad_y.at<float>(u,v));
                }
            }

            cout << "xxx: " << sum1/sum2 << endl;
            // cout << "atan: " << sum1/sum2 << endl;

            Vx.at<float>(i,j) = sum1;
            Vy.at<float>(i,j) = sum2;

            double calc = 0.0;

            if(sum1 != 0 && sum2 != 0)
            {
                calc = 0.5 * atan(sum1 / sum2);
            }

            theta.at<float>(i,j) = calc;

            // Perform low-pass filtering
            float angle = 2 * calc;
            lowPassX.at<float>(i,j) = cos(angle * pi / 180);
            lowPassY.at<float>(i,j) = sin(angle * pi / 180);

            float sum3 = 0.0;
            float sum4 = 0.0;

            for(int u = -lowPassSize / 2; u < lowPassSize / 2; u++)
            {
               for(int v = -lowPassSize / 2; v < lowPassSize / 2; v++)
               {
                  sum3 += inputImage.at<float>(u,v) * lowPassX.at<float>(i - u*lowPassSize, j - v * lowPassSize);
                  sum4 += inputImage.at<float>(u, v) * lowPassY.at<float>(i - u*lowPassSize, j - v * lowPassSize);
               }
            }
        lowPassX2.at<float>(i,j) = sum3;
        lowPassY2.at<float>(i,j) = sum4;

        float calc2 = 0.0;

        if(sum3 != 0 && sum4 != 0)
        {
           calc2 = 0.5 * atan(lowPassY2.at<float>(i, j) / lowPassX2.at<float>(i, j)) * 180 / pi;
        }
        // orientationMat.at<float>(i,j) = calc2;
        cout << calc << endl;
        int length = 8;
        Point P1(i,j);
        Point P2;
        P2.x =  (int)round(P1.x + length * cos(calc2 * CV_PI / 180.0));
        P2.y =  (int)round(P1.y + length * sin(calc2 * CV_PI / 180.0));


        line(orientationMat, P1, P2, CV_RGB(255, 255, 255));
        // imshow("Block", orie);
        }

    }
return orientationMat;

}

int calculate_orientation(Mat &block, Mat &orie) {
    int scale = 1;
    int delta = 0;
    int ddepth = CV_32F;
    Mat grad_x, grad_y;
    Mat abs_grad_x, abs_grad_y;
    Sobel( block, grad_x, ddepth, 1, 0, 3, scale, delta, BORDER_DEFAULT );
    convertScaleAbs( grad_x, abs_grad_x );

    Sobel( block, grad_y, ddepth, 0, 1, 3, scale, delta, BORDER_DEFAULT );
    convertScaleAbs( grad_y, abs_grad_y );

    Mat blocka;

    // if (block_num) {cout << abs_grad_y << endl; cout << abs_grad_x << endl;}
    // cout << block_num << endl;
    imshow("grad_y", abs_grad_y);

    // addWeighted( abs_grad_x, 0.5, abs_grad_y, 0.5, 0, block );

    // imshow("Sobel", blocka);

    // cout << block.rows << " " << block.cols << endl;
    for (int x=block_len/2 ; x < block.rows - block_len/2 ; x += block_len) {
        for (int y=block_len/2 ; y < block.cols - block_len/2 ; y += block_len) {
            float a = 0.0f, b = 0.0f;
            for (int i=x-block_len/2 ; i<x + block_len/2 ; i++) {
                for (int j= y - block_len/2 ; j<y + block_len/2 ; j++) {
                    // cout << float(abs_grad_x.at<uchar>(i,j)) << endl;
                    a += 2*((float)abs_grad_x.at<uchar>(i,j))*((float)abs_grad_y.at<uchar>(i,j));
                    b += ((float)abs_grad_x.at<uchar>(i,j))*((float)abs_grad_x.at<uchar>(i,j))-((float)abs_grad_y.at<uchar>(i,j))*((float)abs_grad_y.at<uchar>(i,j));
                }
             }

            // cout << a << endl;
            // cout << b << endl;
            

            // cout << a/b << endl;
            float angle;
            if (a == 0) {
                angle = 0.0f;
            } else if (b == 0) {
                angle = 90.0f;
            } else {
                angle = 0.5 * atan(a/b) * 180/CV_PI;
                // angle = 15.0f;
            }

            //wieclaw addition (?)
            int k = 0;
            if ((angle < 0 && a<0) || (angle >=0 && a > 0)) {
                k = 90;
            } else if (angle < 0 && a >= 0) {
                k = 180;
            } else if (angle >= 0 && a<= 0) {
                k = 0;
            }

            // angle += k;


            // cout << "Block " << block_num << " : " << angle << endl;

            int length = 8;
            Point P1(x,y);
            if (angle < 0) {
                // angle *= -1;
                // P1.x = 0;
                // P1.y = 16;
            } else {
                // angle *= -1;
                // P1.x = 0;
                // P1.y = 16;
            }
            Point P2;

            // Mat orientation = Mat::zeros(abs_grad_x.rows, abs_grad_y.cols, CV_32F); //to store the gradients grad_x.convertTo(grad_x,CV_32F);

            // grad_x.convertTo(grad_x,CV_32F);

            // grad_y.convertTo(grad_y,CV_32F);

            // phase(grad_x, grad_y, orientation, true);

            // angle = orientation.at<float>(8,8);

            P2.x =  (int)round(P1.x + length * cos(angle * CV_PI / 180.0));
            P2.y =  (int)round(P1.y + length * sin(angle * CV_PI / 180.0));


            line(orie, P1, P2, CV_RGB(255, 255, 255));
            imshow("Block", orie);
            // waitKey(0);
            // block.setTo(Scalar(0,0,0));
            // line(block, P1, P2, CV_RGB(255, 255, 255));
        }
    }
    
    
}

void normalize_image(Mat &image) {
    // float calc_mean;
    // int sum = 0;
    // for (int i=0 ; i<image.rows ; i++) {
    //  for (int j=0 ; j<image.cols ; j++) {
    //      sum += image.at<uchar>(i,j);
    //  }
    // }
    // calc_mean = 1.0f*sum/(image.rows*image.cols);
    // cout << "Calc mean " << calc_mean << endl;

    Scalar mean, stddev;
    meanStdDev(image, mean, stddev);
    cout << "Mean " << mean.val[0] << endl;
    cout << "Stddev " << stddev.val[0] << endl;


    for (int i=0 ; i<image.rows ; i++) {
        for (int j=0 ; j<image.cols ; j++) {
            // if (i>70 && i<90 && j>70 && j<90)cout << i << " " << j << " Before " << (int)image.at<uchar>(i,j);// << endl;
            // if (i < 100 && j<100) cout << i << " " << j << " Before " << (int)image.at<uchar>(i,j);
            int addition = sqrt(100*pow((image.at<uchar>(i,j)-mean.val[0]), 2)/stddev.val[0]);
            if (image.at<uchar>(i,j) > mean.val[0]) {
                if (addition + 100 > 255)
                    image.at<uchar>(i,j) = 255;
                else 
                    image.at<uchar>(i,j) = 100 + addition;
            } else {
                if (addition > 100)
                    image.at<uchar>(i,j) = 0;
                else
                    image.at<uchar>(i,j) = 100 - sqrt(100*pow((image.at<uchar>(i,j)-mean.val[0]), 2)/stddev.val[0]);
            }
            // if (i>70 && i<90 && j>70 && j<90) cout << " After " << (int)image.at<uchar>(i,j) << endl;
            // if ((int)image.at<uchar>(i,j) > 200) 
            //  cout <<  i << " " << j <<" After " << (int)image.at<uchar>(i,j) << endl;
        }
    }
    imshow("Normim", image);
    waitKey(0);
}

int main(int argc, char** argv) {
    //Read image
    String img = argv[1];
    Mat image = imread(img, 0);

    //Check for failure
    if (image.empty()) {
        cout << "Could not open or find the image" << endl;
        return -1;
    }

    Mat normim = image.clone();

    Scalar mean, stddev;
    meanStdDev(image, mean, stddev);
    cout << "Mean " << mean.val[0] << endl;
    cout << "Stddev " << stddev.val[0] << endl;

    // subtract(image, 1.0f*mean.val[0], normim);
    // normim = normim*1.0f/stddev.val[0];
    
    imshow("Fingerprint", image);
    waitKey(0);
    normalize_image(normim);

    meanStdDev(normim, mean, stddev);
    cout << "Mean " << mean.val[0] << endl;
    cout << "Stddev " << stddev.val[0] << endl;

    std::vector<cv::Mat> blocks;

    subdivide(normim, 256/block_len, 256/block_width, blocks);

    for (int i=0 ; i<256/block_len*256/block_width ; i++) {
        mask[i] = segment(blocks[i], 20);
    }

    Mat segmented = image.clone();
    Mat orie = image.clone()  ;
    orie.setTo(Scalar(0,0,0));

    // for (int i=0; i<256/block_len*256/block_width; i++) {
    //     if (mask[i]) calculate_orientation(blocks[i], i);
    // }

    // calculate_orientation(normim, orie);
    orie = orientation(normim);

    imshow("Image", segmented);
    imshow("Orientation", orie);
    // cout << orie << endl;
    waitKey(0);
    return 0;
}

//g++ -ggdb extract_test.cpp -o extract_test `pkg-config --cflags --libs opencv`