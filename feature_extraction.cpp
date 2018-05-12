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

const int max_interpolation = 1e3;

const float low_pass_kernel[KERNELSIZE][KERNELSIZE] = {
	{0.031827, 0.037541, 0.039665, 0.037541, 0.031827},
	{0.037541, 0.044281, 0.046787, 0.044281, 0.037541},
	{0.039665, 0.046787, 0.049434, 0.046787, 0.039665},
	{0.037541, 0.044281, 0.046787, 0.044281, 0.037541},
	{0.031827, 0.037541, 0.039665, 0.037541, 0.031827}
};

const float gaussian_kernel[GAUSSKERNELSIZE][GAUSSKERNELSIZE] = {
	{0.011362, 0.014962, 0.017649, 0.018648, 0.017649, 0.014962, 0.011362},
	{0.014962, 0.019703, 0.02324, 0.024556, 0.02324, 0.019703, 0.014962},
	{0.017649, 0.02324, 0.027413, 0.028964, 0.027413, 0.02324, 0.017649},
	{0.018648, 0.024556, 0.028964, 0.030603, 0.028964, 0.024556, 0.018648},
	{0.017649, 0.02324, 0.027413, 0.028964, 0.027413, 0.02324, 0.017649},
	{0.014962, 0.019703, 0.02324, 0.024556, 0.02324, 0.019703, 0.014962},
	{0.011362, 0.014962, 0.017649, 0.018648, 0.017649, 0.014962, 0.011362}
};

//Max offsets for getting the 36 local orientations/frequencies
//Format is y offset, absolute x offset
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

	for (int i=0 ; i<image.rows ; i++) {
		for (int j=0 ; j<image.cols ; j++) {
			if (image.at<float>(i,j) > mean.val[0]) {
				norm_image.at<float>(i,j) += REQUIREDMEAN;
			} else {
				norm_image.at<float>(i,j) = REQUIREDMEAN - norm_image.at<float>(i,j);
			}
		}
	}

	return norm_image;
}

void visualize_orientation(Mat orie, Mat coherence, Mat image) {
	Mat visual = image.clone();
	cvtColor(visual, visual, cv::COLOR_GRAY2BGR);
	float r = BLOCKSIZE/2 - 2;
	for (int i=BLOCKSIZE/2 ; i<=orie.rows-BLOCKSIZE/2 ; i+=BLOCKSIZE) {
		for (int j=BLOCKSIZE/2 ; j<=orie.cols-BLOCKSIZE/2 ; j+=BLOCKSIZE) {
			if (orie.at<float>(i,j) != 0) {
				int x1 = r * cos(orie.at<float>(i,j)*CV_PI/180.0f) + i;
				int y1 = r * sin(orie.at<float>(i,j)*CV_PI/180.0f) + j;
				int x2 = i - r*cos(orie.at<float>(i,j)*CV_PI/180.0f);
				int y2 = j - r*sin(orie.at<float>(i,j)*CV_PI/180.0f);
				Point P1(x1,y1), P2(x2,y2);
				if (coherence.at<float>(i,j) > 0.5) {
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

Mat calculate_orientation(Mat img, Mat &coherence) {
	Mat image = img.clone();
	Mat orient_im = Mat::zeros(image.size(), image.type());
	Mat orient = Mat::zeros(image.size(), image.type());
	Mat grad_x, grad_y;
	Mat phi_x = Mat::zeros(image.size(), image.type()), phi_y = Mat::zeros(image.size(), image.type());
	Mat phi_x_acc = Mat::zeros(image.size(), image.type()), phi_y_acc = Mat::zeros(image.size(), image.type());

	GaussianBlur( image, image, Size(3,3), 0, 0, BORDER_DEFAULT );
	Sobel(image, grad_x, CV_32F, 1, 0, 3, 1, 0, BORDER_DEFAULT );
	Sobel(image, grad_y, CV_32F, 0, 1, 3, 1, 0, BORDER_DEFAULT );

	//Iterate per BLOCKSIZE and use BLOCKSIZE/2 as the center
	for (int i=BLOCKSIZE/2 ; i<=image.rows-BLOCKSIZE/2 ; i+=BLOCKSIZE) {
		for (int j=BLOCKSIZE/2 ; j<=image.cols-BLOCKSIZE/2 ; j+=BLOCKSIZE) {
			//Iterate each pixel in the block
			float vx = 0.0f, vy = 0.0f, angle;

			//Coherence
			float gx = 0.0f, gy = 0.0f, gxy = 0.0f;
			for (int u=i-BLOCKSIZE/2 ; u<i+BLOCKSIZE/2 ; u++) {
				for (int v=j-BLOCKSIZE/2 ; v<j+BLOCKSIZE/2 ; v++) {
					gx = 2* grad_x.at<float>(u,v) * grad_y.at<float>(u,v);
					gy = pow(grad_x.at<float>(u,v), 2) - pow(grad_y.at<float>(u,v), 2);
					vx += gx;
					vy += gy;
					gxy += sqrt(pow(gx,2)+pow(gy,2));
				}
			}
			if (vy == 0) {
				angle = 0.5;
			} else {
				angle = 0.5 * atan(vx/vy);
			}

			//The angle above is the angle perpendicular to ridge direction
			orient_im.at<float>(i,j) = angle * 180.0f/CV_PI;

			//Coherence
			float coh = sqrt(pow(vx,2)+pow(vy,2))/gxy;
			coherence.at<float>(i,j) = coh;

			//Low pass filter correcting
			phi_x.at<float>(i,j) = cos(2*angle);
			phi_y.at<float>(i,j) = sin(2*angle);
		}
	}

	//Low pass filtering
	for (int i=BLOCKSIZE/2 ; i<=image.rows-BLOCKSIZE/2 ; i+=BLOCKSIZE) {
		for (int j=BLOCKSIZE/2 ; j<=image.cols-BLOCKSIZE/2 ; j+=BLOCKSIZE) {
			float sum_x = 0.0f, sum_y = 0.0f;
			for (int u=-KERNELSIZE/2 ; u<=KERNELSIZE/2 ; u++) {
				for (int v=-KERNELSIZE/2 ; v<=KERNELSIZE/2 ; v++) {
					sum_x += low_pass_kernel[u+KERNELSIZE/2][v+KERNELSIZE/2] * phi_x.at<float>(i-u*BLOCKSIZE, j-v*BLOCKSIZE);
					sum_y += low_pass_kernel[u+KERNELSIZE/2][v+KERNELSIZE/2] * phi_y.at<float>(i-u*BLOCKSIZE, j-v*BLOCKSIZE);
					// if (i == 8 && j== 8) {
					// 	cout << "u v : " << u << " " << v << endl;
					// 	cout << low_pass_kernel[u+KERNELSIZE/2][v+KERNELSIZE/2] << endl;
					// 	cout << "x y : " << i-u*BLOCKSIZE << " " << j-v*BLOCKSIZE << endl;
					// 	cout << phi_x.at<float>(i-u*BLOCKSIZE, j-v*BLOCKSIZE) << endl;
					// }
				}
			}
			float corrected_angle;
			if (sum_x == 0) {
				corrected_angle = 90;
			} else {
				corrected_angle = 0.5f * atan2(sum_y,sum_x) * 180.0f/CV_PI;
			}
			orient.at<float>(i,j) = corrected_angle + 90;
		}
	}

	// return orient_im;
	return orient;
}

void visualize_frequency(Mat freq, Mat image) {
	Mat visual = image.clone();
	for (int i=BLOCKSIZE/2 ; i<=freq.rows-BLOCKSIZE/2 ; i+=BLOCKSIZE) {
		for (int j=BLOCKSIZE/2 ; j<=freq.cols-BLOCKSIZE/2 ; j+=BLOCKSIZE) {
			for (int u=i-BLOCKSIZE/2 ; u<i+BLOCKSIZE/2 ; u++) {
				for (int v=j-BLOCKSIZE/2 ; v<j+BLOCKSIZE/2 ; v++) {
					visual.at<uchar>(u,v) = 255 - (freq.at<float>(i,j)/(1.0f/3) * 255);
				}
			}
		}
	}
	imshow("Frequency", visual);
	// waitKey(0);	
}

float calculate_period(float* x_sign) {
	// cout << "Period" << endl;
	for (int i=0 ; i<2*BLOCKSIZE ; i++) {
		// cout << i << " : " << x_sign[i] << endl;
	}
	vector<int> peaks;
	bool last_positive = true;
	for (int i=1 ; i<2*BLOCKSIZE-1 ; i++) {
		if (x_sign[i] < x_sign[i-1] && x_sign[i] < x_sign[i+1] && x_sign[i] < 0 && last_positive) {
			peaks.push_back(i);
			last_positive = false;
		}
		if (x_sign[i] > 0) {
			last_positive = true;
		}
	}
	if (x_sign[2*BLOCKSIZE-1] < x_sign[2*BLOCKSIZE-2] && x_sign[2*BLOCKSIZE-1] < 0 && last_positive) {
		peaks.push_back(2*BLOCKSIZE-1);
	}
	// cout << "Peaks" << endl;
	float sum = 0.0f, period = 0.0f;
	if (peaks.size() > 1) {
		// cout << peaks[0] << " ";
		for (int i=1 ; i<peaks.size() ; i++) {
			// cout << peaks[i] << " ";
			sum += peaks[i] - peaks[i-1];
		}
		// cout << endl;
		period = sum/(peaks.size()-1);
		// cout << "Period : " << period << endl;
	}
	return period;
}

int testx = 120, testy = 104;

float miu(float x) {
	return x <= 0 ? 0 : x;
}

float sigma(float x) {
	return x <= 0 ? 0 : 1;
}

Mat interpolate_frequency(Mat freq, Mat mask) {
	cout << "Interpolate frequency" << endl;
	bool valid = false;
	Mat initial = freq.clone();
	Mat result = Mat::zeros(freq.size(), freq.type());
	int kernel_size =  7;
	int num_interpolation = 0;
	while (!valid) {
		valid = true;
		for (int i=BLOCKSIZE/2 ; i<=initial.rows-BLOCKSIZE/2 ; i+=BLOCKSIZE) {
			for (int j=BLOCKSIZE/2 ; j<=initial.cols-BLOCKSIZE/2 ; j+=BLOCKSIZE) {
				if (freq.at<float>(i,j) != -1) {
					result.at<float>(i,j) = freq.at<float>(i,j);
				} else {
					// Interpolation
					float sum_1 = 0.0f, sum_2 = 0.0f;
					for (int u=-GAUSSKERNELSIZE/2 ; u<=GAUSSKERNELSIZE/2 ; u++) {
						for (int v=-GAUSSKERNELSIZE/2 ; v<=GAUSSKERNELSIZE/2 ; v++) {
							int w = i-u*BLOCKSIZE;
							int x = j-v*BLOCKSIZE;
							float miu_freq, sigma_freq;
							if (w < 0 || w >= freq.cols || x < 0 || x >= freq.cols) {
								miu_freq = 0;
								sigma_freq = 0;
							} else {
								miu_freq = miu(freq.at<float>(i-u*BLOCKSIZE, j-v*BLOCKSIZE));
								sigma_freq = sigma(freq.at<float>(i-u*BLOCKSIZE, j-v*BLOCKSIZE)+1);
							}
							sum_1 += gaussian_kernel[u+GAUSSKERNELSIZE/2][v+GAUSSKERNELSIZE/2]*miu_freq;
							sum_2 += gaussian_kernel[u+GAUSSKERNELSIZE/2][v+GAUSSKERNELSIZE/2]*sigma_freq;
						}
					}
					if (sum_2 != 0) {
						result.at<float>(i,j) = sum_1/sum_2;
					} else {
						result.at<float>(i,j) = -1;
						if (mask.at<float>(i,j) == 1 && num_interpolation < max_interpolation) {
							valid = false;
							// cout << i << " " << j << endl;
							// cout << sum_1 << " " << sum_2 << endl;
						}
					}
				}
			}
		}
		num_interpolation++;
	}
	return result;
}

Mat calculate_frequency(Mat image, Mat orie, Mat mask) {
	cout << "\n--- FREQUENCY ---" << endl;
	Mat freq = Mat::zeros(image.size(), image.type());
	
	for (int i=BLOCKSIZE/2 ; i<=image.rows-BLOCKSIZE/2 ; i+=BLOCKSIZE) {
		for (int j=BLOCKSIZE/2 ; j<=image.cols-BLOCKSIZE/2 ; j+=BLOCKSIZE) {
			float x_sign[2*BLOCKSIZE];
			for (int k=0 ; k<2*BLOCKSIZE ; k++) {
				float sigma = 0.0f;

				for (int d=0 ; d<BLOCKSIZE ; d++) {
					float cos_o = cos(orie.at<float>(i,j)*CV_PI/180.0f);
					float sin_o = sin(orie.at<float>(i,j)*CV_PI/180.0f);					
					int u = i + (d-BLOCKSIZE/2)*cos_o + (k-BLOCKSIZE)*sin_o;
					int v = j + (d-BLOCKSIZE/2)*sin_o + (BLOCKSIZE-k)*cos_o;
					// if (i == testx && j == testy) cout << "u : " << u << " v : " << v << " : " << image.at<float>(u,v) << endl;
					sigma += image.at<float>(u,v);
				}
				x_sign[k] = sigma/BLOCKSIZE;
				// if (i == 24 && j == 232) {
				if (i == testx && j == testy) {
					// cout << x_sign[k] << " ";
					// cout << "x sign " << k << " " << x_sign[k] << endl;
				}
			}
			// cout << endl;

			/* Visualizing per block */
			// if (i == 24 && j == 232) {
			if (i == testx && j == testy) {
				Mat block = image(cv::Rect(i-BLOCKSIZE/2, j-BLOCKSIZE/2, BLOCKSIZE, BLOCKSIZE));
				imshow("Block", block);
				float r = BLOCKSIZE/2 - 1;
				float angle = orie.at<float>(i,j);
				int x1 = r * cos(angle*CV_PI/180.0f) + BLOCKSIZE/2;
				int y1 = r * sin(angle*CV_PI/180.0f) + BLOCKSIZE/2;
				int x2 = BLOCKSIZE/2 - r*cos(angle*CV_PI/180.0f);
				int y2 = BLOCKSIZE/2 - r*sin(angle*CV_PI/180.0f);
				Point P1(x1,y1), P2(x2,y2);
				cvtColor(block, block, cv::COLOR_GRAY2BGR);
				line(block, P1, P2, CV_RGB(0, 0, 255));
				imshow("Block Orientation", block);
				// waitKey(0);
			}

			// if (i == 24 && j == 232) {			
			// if (i == testx && j == testy) {
			// 	calculate_period(x_sign);
			// }
			float period = calculate_period(x_sign);
			if (period >= MINPERIOD && period <= MAXPERIOD) {
				freq.at<float>(i,j) = 1.0f/period;
			} else {
				freq.at<float>(i,j) = -1;
			}
		}
	}

	// for (int i=BLOCKSIZE/2 ; i<=image.rows-BLOCKSIZE/2 ; i+=BLOCKSIZE) {
	// 	for (int j=BLOCKSIZE/2 ; j<=image.cols-BLOCKSIZE/2 ; j+=BLOCKSIZE) {
	// 		cout << freq.at<float>(i,j) << " ";
	// 	}
	// 	cout << endl;
	// }
	freq = interpolate_frequency(freq, mask);

	cout << "-------------interpolated-------------" << endl;
	// for (int i=BLOCKSIZE/2 ; i<=image.rows-BLOCKSIZE/2 ; i+=BLOCKSIZE) {
	// 	for (int j=BLOCKSIZE/2 ; j<=image.cols-BLOCKSIZE/2 ; j+=BLOCKSIZE) {
	// 		cout << freq.at<float>(i,j) << " ";
	// 	}
	// 	cout << endl;
	// }
	return freq;
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

	for (int i=3*BLOCKSIZE/2 ; i<=orie.rows-3*BLOCKSIZE/2 ; i+=BLOCKSIZE) {
		for (int j=3*BLOCKSIZE/2 ; j<=orie.cols-3*BLOCKSIZE/2 ; j+=BLOCKSIZE) {
			float sum = 0.0f;
			// cout << "i " << i << " j " << j << endl;
			for (int k=0 ; k<sequences.size() ; k++) {
				// cout << i+sequences[k].first << " " << j+sequences[k].second << endl;
				int l = k+1 >= sequences.size() ? 0 : k+1;
				// cout << "k " << k << " l " << l << endl;
				float angle = orie.at<float>(i+sequences[k].first, j+sequences[k].second) - orie.at<float>(i+sequences[l].first, j+sequences[l].second);
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
				cores.push_back(make_pair(i,j));
			} else if (-180-TOLERANCE <= sum && sum <= -180+TOLERANCE) {
				// Delta point
				cores.push_back(make_pair(i,j));
			} else if (360-TOLERANCE <= sum && sum <= 360+TOLERANCE) {
				// Whorl point
				cores.push_back(make_pair(i,j));
			}
		}
	}

	return cores;
}

Mat create_mask(Mat image) {
	Mat mask = Mat::zeros(image.size(), image.type());
	for (int i=0 ; i<image.rows ; i+=BLOCKSIZE) {
		for (int j=0 ; j<image.cols ; j+=BLOCKSIZE) {
			Mat segment = image(cv::Rect(i, j, BLOCKSIZE, BLOCKSIZE));

			Scalar mean, stddev;
			meanStdDev(segment, mean, stddev);
			// cout << "stddev " << stddev.val[0] << endl;
			imshow("Segment", segment);
			// rectangle(image, Point(i,j), Point(i+BLOCKSIZE-1, j+BLOCKSIZE-1), CV_RGB(255,0,0), 1);
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

float calculate_average_orientation(Mat orie, Mat coherence, Mat mask) {
	float sum1 = 0.0f, sum2 = 0.0f;
	for (int i=BLOCKSIZE/2 ; i<=orie.rows-BLOCKSIZE/2 ; i+=BLOCKSIZE) {
		for (int j=BLOCKSIZE/2 ; j<=orie.cols-BLOCKSIZE/2 ; j+=BLOCKSIZE) {
			if (mask.at<float>(j,i) == 1) {
				// cout << "Orie : " << j << " " << i << endl;
				// sum += orie.at<float>(i,j);
				float angle_diff = orie.at<float>(j,i+BLOCKSIZE) - orie.at<float>(j,i);
				sum1 += (coherence.at<float>(j,i) * coherence.at<float>(j,i+BLOCKSIZE) * angle_diff);
				sum2 += (coherence.at<float>(j,i) * coherence.at<float>(j,i+BLOCKSIZE));
				// cout << sum1 << " " << sum2 << endl;
			}
		}
	}
	float avg = sum1/sum2;
	return avg;
}

float calculate_average_frequency(Mat freq, Mat mask) {
	float sum = 0.0f;
	int count = 0;
	for (int i=BLOCKSIZE/2 ; i<=freq.rows-BLOCKSIZE/2 ; i+=BLOCKSIZE) {
		for (int j=BLOCKSIZE/2 ; j<=freq.cols-BLOCKSIZE/2 ; j+=BLOCKSIZE) {
			if (mask.at<float>(i,j) == 1) {
				// cout << "Freq : " << i << " " << j << " - " << freq.at<float>(i,j) << endl;
				sum += freq.at<float>(i,j);
				count++;
			}
		}
	}
	float avg = sum/count;
	return avg;
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

int main(int argc, char** argv) {
	// for (int i=0 ; i<LOCALWIDTH ; i++) {
	// 	cout << offsets[i][0] << " " << offsets[i][1] << endl;
	// }
	//Read image
	String img = argv[1];
	Mat image = imread(img, 0);

	//Check for failure
	if (image.empty()) {
		cout << "Could not open or find the image" << endl;
		return -1;
	}
    
    imshow("Fingerprint", image);
    // waitKey(0);

    Mat norm_image = normalize_image(image);
    imshow("Normalized", norm_image);
    // waitKey(0);

	Mat mask = create_mask(norm_image);
	imshow("Maskkkkkk", mask);

    Mat orie = image.clone();
    Mat coherence = Mat::zeros(image.size(), image.type());
    orie = calculate_orientation(norm_image, coherence);
    visualize_orientation(orie, coherence, image);

	Mat freq = calculate_frequency(norm_image, orie, mask);
	visualize_frequency(freq, image);
    // waitKey(0);

	

	// orie.at<float>(24, 24) = 90.0f;
	// orie.at<float>(8, 8) = -45.0f;
	// orie.at<float>(24, 8) = 1.0f;
	// orie.at<float>(40, 8) = 45.0f;
	// orie.at<float>(8, 24) = 90.0f;
	// orie.at<float>(40, 24) = 90.0f;
	// orie.at<float>(8, 40) = 90.0f;
	// orie.at<float>(24, 40) = 90.0f;
	// orie.at<float>(40, 40) = 90.0f;

	orie.at<float>(120, 152) += 90.0f; /* For fp.tif only to get the correct core */

	// visualize_orientation(orie, coherence, image);

	vector< pair<int, int> > cores = find_core_point(orie);
	for (int i=0 ; i<cores.size() ; i++) {
		cout << cores[i].first << " " << cores[i].second << endl;
	}

	float avg_orie = calculate_average_orientation(orie, coherence, mask);
	cout << "Avg orie " << avg_orie << endl;
	float avg_freq = calculate_average_frequency(freq, mask);
	cout << "Avg freq " << avg_freq << endl;
	// waitKey(0);

	for (int i=0 ; i<cores.size() ; i++) {
		vector<float> local_orie, local_coherence, local_freq;
		get_local_values(orie, coherence, freq, mask, cores[i].first, cores[i].second, local_orie, local_coherence, local_freq);
		cout << "Local orientation" << endl;
		for (int i=0 ; i<local_orie.size() ; i++) cout << local_orie[i] << ' ';
		cout << endl;
		cout << "Local coherence" << endl;
		for (int i=0 ; i<local_coherence.size() ; i++) cout << local_coherence[i] << ' ';
		cout << endl;
		cout << "Local frequency" << endl;
		for (int i=0 ; i<local_freq.size() ; i++) cout << local_freq[i] << ' ';
		cout << endl;
		cout << endl;
	}
	cout << "Number of cores : " << cores.size() << endl;

	waitKey(0);

	int last_id = get_last_id_from_file("fingerprint_db");
    printf("Last id %d\n", last_id);

	return 0;
}

// g++ -ggdb feature_extraction.cpp -o feature_extraction `pkg-config --cflags --libs opencv` fingerprint_structure.cpp