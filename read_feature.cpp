#include <opencv2/opencv.hpp>
#include <iostream>
#include <vector>
#include "fingerprint_structure.hpp"
using namespace cv;
using namespace std;

#define THRESHOLD 40
#define BLOCKSIZE 16
#define REQUIREDMEAN 100
#define REQUIREDVAR 100
#define MINPERIOD 3
#define MAXPERIOD 25
#define KERNELSIZE 5
#define GAUSSKERNELSIZE 7
#define TOLERANCE 50
#define LOCALWIDTH 6

const int offsets[LOCALWIDTH][2] = {{-1, 3}, {0,4}, {1,3}, {2,3}, {3,2}, {4,0}};

Mat normalize_image(Mat image) {
	Mat norm_image = image.clone();
	norm_image.convertTo(norm_image, CV_32FC1);

	Scalar mean, stddev;
    meanStdDev(norm_image, mean, stddev);
    cout << "Mean " << mean.val[0] << endl;
    cout << "Stddev " << stddev.val[0] << endl;

	norm_image = norm_image - mean.val[0];
	pow(norm_image, 2, norm_image);
	norm_image = norm_image * REQUIREDVAR / stddev.val[0];
	sqrt(norm_image, norm_image);

	for (int i=0 ; i<image.cols ; i++) {
		for (int j=0 ; j<image.rows ; j++) {
			if (image.at<float>(j,i) > mean.val[0]) {
				norm_image.at<float>(j,i) += REQUIREDMEAN;
			} else {
				norm_image.at<float>(j,i) = REQUIREDMEAN - norm_image.at<float>(j,i);
			}
		}
	}

	return norm_image;
}

Mat create_mask(Mat image) {
	Mat mask = Mat::zeros(image.size(), image.type());
	for (int i=0 ; i<image.cols ; i+=BLOCKSIZE) {
		for (int j=0 ; j<image.rows ; j+=BLOCKSIZE) {
			Mat segment = image(cv::Rect(i, j, BLOCKSIZE, BLOCKSIZE));

			Scalar mean, stddev;
			meanStdDev(segment, mean, stddev);
			// cout << "stddev " << stddev.val[0] << endl;
			// imshow("Segment", segment);
			// rectangle(image, Point(j,i), Point(i+BLOCKSIZE-1, j+BLOCKSIZE-1), CV_RGB(255,0,0), 1);
			// imshow("Image segment", image);
			// waitKey(0);


			// If stddev is more than threshold then it is fingerprint area
			if (stddev.val[0] > THRESHOLD) {
				// cout << "Masked : " << i << " " << j << endl;
				for (int u=0 ; u<BLOCKSIZE ; u++) {
					for (int v=0 ; v<BLOCKSIZE ; v++) {
						mask.at<float>(j+u,i+v) = 1;
					}
				}
			} else {
				for (int u=0 ; u<BLOCKSIZE ; u++) {
					for (int v=0 ; v<BLOCKSIZE ; v++) {
						image.at<float>(j+u,i+v) = 0;
					}
				}
			}
		}
	}
	// imshow("Mask", mask);
	// imshow("Masked image", image);
	// waitKey(0);
	return mask;
}

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

vector< pair<int,int> > find_core_point(Mat orie) {
	cout << "CORE" << endl;
	// Initiate location vectors
	vector< pair<int, int> > sequences;
	sequences.push_back(make_pair(-BLOCKSIZE, BLOCKSIZE));
	sequences.push_back(make_pair(-BLOCKSIZE, 0));
	sequences.push_back(make_pair(-BLOCKSIZE, -BLOCKSIZE));
	sequences.push_back(make_pair(0, -BLOCKSIZE));
	sequences.push_back(make_pair(BLOCKSIZE, -BLOCKSIZE));
	sequences.push_back(make_pair(BLOCKSIZE, 0));
	sequences.push_back(make_pair(BLOCKSIZE, BLOCKSIZE));
	sequences.push_back(make_pair(0, BLOCKSIZE));

	vector< pair<int,int> > cores;

	for (int i=3*BLOCKSIZE/2 ; i<=orie.cols-3*BLOCKSIZE/2 ; i+=BLOCKSIZE) {
		for (int j=3*BLOCKSIZE/2 ; j<=orie.rows-3*BLOCKSIZE/2 ; j+=BLOCKSIZE) {
			float sum = 0.0f;
			// cout << "i " << i << " j " << j << endl;
			for (int k=0 ; k<sequences.size() ; k++) {
				// cout << i+sequences[k].first << " " << j+sequences[k].second << endl;
				int l = k+1 >= sequences.size() ? 0 : k+1;
				// cout << "k " << k << " l " << l << endl;
				float angle = orie.at<float>(j+sequences[k].first, i+sequences[k].second) - orie.at<float>(j+sequences[l].first, i+sequences[l].second);
				// cout << angle << endl;
				if (angle > 90) {
					angle -= 180;
				} else if (angle < -90) {
					angle += 180;
				}
				sum += angle;
			}
			// cout << sum << endl;
			if (180-TOLERANCE <= sum && sum <= 180+TOLERANCE) {
				// Loop point
				cores.push_back(make_pair(j,i));
			} else if (-180-TOLERANCE <= sum && sum <= -180+TOLERANCE) {
				// Delta point
				cores.push_back(make_pair(j,i));
			} else if (360-TOLERANCE <= sum && sum <= 360+TOLERANCE) {
				// Whorl point
				cores.push_back(make_pair(j,i));
			}
		}
	}

	return cores;
}

void get_local_values(Mat orie, Mat coherence, Mat freq, Mat mask, int core_i, int core_j, vector<float> &local_orie, vector<float> &local_coherence, vector<float> &local_freq) {
	cout << "CORE " << core_i << " " << core_j << endl;
	for (int i=0 ; i<LOCALWIDTH ; i++) {
		int y_offset = offsets[i][0];
		int x_offset = offsets[i][1];
		for (int j=-x_offset ; j<=x_offset ; j++) {
			int u = core_j+y_offset*BLOCKSIZE;
			int v = core_i+j*BLOCKSIZE;
			// cout << u << " " << v << endl;
			// mask.at<float>(u, v) = 0;
			if (mask.at<float>(u,v) != 0) {
				local_orie.push_back(orie.at<float>(u,v));
				local_coherence.push_back(coherence.at<float>(u,v));
				local_freq.push_back(freq.at<float>(u,v));
			} else {
				local_orie.push_back(0);
				local_coherence.push_back(0);
				local_freq.push_back(0);
			}
		}
	}
	// imshow("Locals", mask);
	// waitKey(0);
}

void visualize_core(Mat orie, Mat coherence, Mat image, const vector< pair<int, int> > &cores, Mat mask) {
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
    int j = 0;
    for (int i=0 ; i<cores.size() && j<5 ; i++) {
        // visual.at<Vec3b>(cores[i].first, cores[i].second) = Vec3b(0, 255, 0);
        vector<float> local_orie, local_coherence, local_freq;
        get_local_values(orie, coherence, orie, mask, cores[i].first, cores[i].second, local_orie, local_coherence, local_freq);
        int notzero = 0;
        for (int k=0 ; k<36 ; k++) {
            if (local_orie[k] != 0) {
                notzero += 1;
                if (notzero >= 18) break;
            }
        }
        if (notzero >= 18) {
            circle(visual, Point(cores[i].second, cores[i].first), 4, Vec3b(0, 255, 0), -1);
            j++;
        }
    }
	imshow("Cores", visual);
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

    Mat norm_image = normalize_image(img);

	Mat mask = create_mask(norm_image);
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
    vector< pair<int, int> > cores = find_core_point(orie);

    // visualize_orientation(orie, cohe, img);
    // waitKey(0);
    // visualize_frequency(freq, img);
    // waitKey(0);
    visualize_core(orie, cohe, img, cores, mask);
    waitKey(0);
    cout << orie.depth() << " " << orie.channels() << endl;
    cout << cohe.depth() << " " << cohe.channels() << endl;
    cout << freq.depth() << " " << freq.channels() << endl;
    return 0;
}

// g++ -ggdb read_feature.cpp -o read_feature `pkg-config --cflags --libs opencv` fingerprint_structure.cpp -std=c++11
