#include <opencv2/opencv.hpp>
#include <iostream>
#include <vector>
#include "fingerprint_structure.hpp"
using namespace cv;
using namespace std;

#define BLOCKSIZE 16

void visualize_orientation(Mat orie, Mat coherence, Mat image) {
	Mat visual = image.clone();
	cvtColor(visual, visual, cv::COLOR_GRAY2BGR);
	float r = BLOCKSIZE/2 - 2;
	for (int i=BLOCKSIZE/2 ; i<=orie.cols-BLOCKSIZE/2 ; i+=BLOCKSIZE) {
		for (int j=BLOCKSIZE/2 ; j<=orie.rows-BLOCKSIZE/2 ; j+=BLOCKSIZE) {
			
			if (orie.at<float>(j,i) != 0) {
				// cout << j << " " << i << " " << orie.at<float>(j,i) << endl; 
				int x1 = r * cos(orie.at<float>(j,i)*CV_PI/180.0f) + i;
				int y1 = r * sin(orie.at<float>(j,i)*CV_PI/180.0f) + j;
				int x2 = i - r*cos(orie.at<float>(j,i)*CV_PI/180.0f);
				int y2 = j - r*sin(orie.at<float>(j,i)*CV_PI/180.0f);
				// cout << j << " " << i << " " << orie.at<float>(j,i) << " " << x1 << " " << y1 << " " << x2 << " " << y2 << " " << endl; 
				Point P1(x1,y1), P2(x2,y2);
				if (coherence.at<float>(j,i) > 0.5) {
					line(visual, P1, P2, CV_RGB(0, 0, 255));
				} else {
					line(visual, P1, P2, CV_RGB(255, 0, 0));
				}
			}
		}
	}
	imshow("Orientation", visual);
	// waitKey(0);	
}

void visualize_frequency(Mat freq, Mat image) {
	Mat visual = image.clone();
	for (int i=BLOCKSIZE/2 ; i<=freq.cols-BLOCKSIZE/2 ; i+=BLOCKSIZE) {
		for (int j=BLOCKSIZE/2 ; j<=freq.rows-BLOCKSIZE/2 ; j+=BLOCKSIZE) {
			// cout << j << " " << i << " " << freq.at<float>(j,i) << endl;
			for (int u=i-BLOCKSIZE/2 ; u<i+BLOCKSIZE/2 ; u++) {
				for (int v=j-BLOCKSIZE/2 ; v<j+BLOCKSIZE/2 ; v++) {
					visual.at<uchar>(v,u) = 255 - (freq.at<float>(j,i)/(1.0f/3) * 255);
				}
			}
		}
	}
	imshow("Frequency", visual);
	// waitKey(0);	
}

int main(int argc, char** argv) {
    if (argc < 2) {
        cerr << "Usage : ./read_feature name\n";
        return 0;
    }

	//Read image
	String name = argv[1];
    string filename = name.substr(0, name.find("."));
    Mat img = imread(name, 0);
    /*Mat orie = imread("fea/"+filename+"_orie.png", 0);
    orie.convertTo(orie,CV_32F);
    Mat cohe = imread("fea/"+filename+"_cohe.png", 0);
    // cohe.convertTo(cohe,CV_32F);
    Mat freq = imread("fea/"+filename+"_freq.png", 0);
    freq.convertTo(freq, CV_32F);
    */
   cout << " a\n";
    FileStorage fs("fea/"+filename+"_fea.yml", FileStorage::READ);
    cout << "b\n";
    Mat orie;
    Mat cohe;
    Mat freq;
    fs["orie"] >> orie;
    fs["cohe"] >> cohe;
    fs["freq"] >> freq;
    cout << "c\n";
    for (int i=0 ; i<orie.rows ; i++) {
        for (int j=0 ; j<orie.cols ; j++) {
            // cout << orie.at<float>(i,j) << " ";
            cout << cohe.at<float>(i,j) << " ";
            // freq.at<float>(i,j) = freq.at<float>(i,j)/100000;
            // cout << freq.at<float>(i,j) << " ";
        }
        cout << endl;
    }

    visualize_orientation(orie, cohe, img);
    waitKey(0);
    visualize_frequency(freq, img);
    waitKey(0);
    cout << orie.depth() << " " << orie.channels() << endl;
    cout << cohe.depth() << " " << cohe.channels() << endl;
    cout << freq.depth() << " " << freq.channels() << endl;
    return 0;
}

// g++ -ggdb read_feature.cpp -o read_feature `pkg-config --cflags --libs opencv` fingerprint_structure.cpp
