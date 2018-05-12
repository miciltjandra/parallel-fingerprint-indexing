#include <stdio.h>
#include <stdarg.h>
#include "opencv2/opencv.hpp"
using namespace std;
using namespace cv;
int main(int argc, char* argv[])
{
    namedWindow("source");
    namedWindow("result");

    namedWindow("ang");

    Mat img=imread("fp.tif",0);
    cv::threshold(img,img,128,255,cv::THRESH_BINARY);
    Mat thinned;

    thinned=img.clone(); // Just clone the input
    //Thinning(img,thinned); // Not actually needed

    cv::GaussianBlur(thinned,thinned,Size(3,3),1.0);
    Mat gx,gy,ang,mag;
    cv::Sobel(thinned,gx,CV_32FC1,1,0);
    cv::Sobel(thinned,gy,CV_32FC1,0,1);
    cv::phase(gx,gy,ang,false);
    cv::magnitude(gx,gy,mag);

    cv::normalize(mag,mag,0,1,cv::NORM_MINMAX);


    Mat angRes=Mat::zeros(img.rows*3,img.cols*3,CV_8UC1);

    for (int i=0;i< img.rows;i+=2)
    {
        for (int j=0;j< img.cols;j+=2)
        {
            int x=j*3;
            int y=i*3;

            float r=5;
            float m=r*(mag.at<float>(i,j));
            float dx=m*r*cos(ang.at<float>(i,j));
            float dy=m*r*sin(ang.at<float>(i,j));

            cv::line(angRes,cv::Point(x,y),cv::Point(x+dx,y+dy),Scalar::all(255),1,CV_AA);
        }
    }
    imshow("ang",angRes);
    imshow("source",img);
    imshow("result",thinned);
    cv::waitKey(0);
}

//g++ -ggdb another_test.cpp -o another_test `pkg-config --cflags --libs opencv`