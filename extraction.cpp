#include <opencv2/opencv.hpp>
#include <iostream>
using namespace cv;
using namespace std;

#define THRESHOLD 10
#define block_len 16
#define block_width 16

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

int calculate_orientation(Mat &block, int block_num) {
	int scale = 1;
	int delta = 0;
	int ddepth = CV_16S;
	Mat grad_x, grad_y;
	Mat abs_grad_x, abs_grad_y;
	Sobel( block, grad_x, ddepth, 1, 0, 3, scale, delta, BORDER_DEFAULT );
	convertScaleAbs( grad_x, abs_grad_x );

	Sobel( block, grad_y, ddepth, 0, 1, 3, scale, delta, BORDER_DEFAULT );
  	convertScaleAbs( grad_y, abs_grad_y );

  	Mat blocka;

  	// if (block_num) {cout << abs_grad_y << endl; cout << abs_grad_x << endl;}
  	cout << block_num << endl;
  	imshow("grad_y", abs_grad_y);

  	// addWeighted( abs_grad_x, 0.5, abs_grad_y, 0.5, 0, block );

  	// imshow("Sobel", blocka);

  	// cout << block.rows << " " << block.cols << endl;
  	float a = 0.0f, b = 0.0f;
  	for (int i=0 ; i<block.rows ; i++) {
  		for (int j=0 ; j<block.cols ; j++) {
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


  	cout << "Block " << block_num << " : " << angle << endl;

	int length = 8;
	Point P1(8,8);
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


	line(block, P1, P2, CV_RGB(0, 255, 0));
	imshow("Block", block);
	// waitKey(0);
	block.setTo(Scalar(0,0,0));
	line(block, P1, P2, CV_RGB(255, 255, 255));
	
}

void normalize_image(Mat &image) {
	// float calc_mean;
	// int sum = 0;
	// for (int i=0 ; i<image.rows ; i++) {
	// 	for (int j=0 ; j<image.cols ; j++) {
	// 		sum += image.at<uchar>(i,j);
	// 	}
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
			// 	cout <<  i << " " << j <<" After " << (int)image.at<uchar>(i,j) << endl;
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

	for (int i=0; i<256/block_len*256/block_width; i++) {
		if (mask[i]) calculate_orientation(blocks[i], i);
	}

	imshow("Image", segmented);
	imshow("Orientation", normim);
	waitKey(0);
	return 0;
}

//g++ -ggdb extraction.cpp -o extraction `pkg-config --cflags --libs opencv`