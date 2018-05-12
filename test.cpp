#include <opencv2/opencv.hpp>
#include <iostream>
using namespace cv;
using namespace std;

int mask[256];

// function for splitting image into multiple blocks. rowDivisor and colDivisor specify the number of blocks in rows and cols respectively
int subdivide(cv::Mat &img, const int rowDivisor, const int colDivisor, std::vector<cv::Mat> &blocks)
{        
    /* Checking if the image was passed correctly */
    if(!img.data || img.empty())
        std::cerr << "Problem Loading Image" << std::endl;

    /* Cloning the image to another for visualization later, if you do not want to visualize the result just comment every line related to visualization */
    cv::Mat maskImg = img.clone();
    /* Checking if the clone image was cloned correctly */
    if(!maskImg.data || maskImg.empty())
        std::cerr << "Problem Loading Image" << std::endl;

    // check if divisors fit to image dimensions
    if(img.cols % colDivisor == 0 && img.rows % rowDivisor == 0)
    {
    	int k =0;
        for(int y = 0; y < img.cols; y += img.cols / colDivisor)
        {
            for(int x = 0; x < img.rows; x += img.rows / rowDivisor)
            {
                blocks.push_back(img(cv::Rect(y, x, (img.cols / colDivisor), (img.rows / rowDivisor))));
                Scalar mean, stddev;
				meanStdDev(blocks[k], mean, stddev);

				if (stddev.val[0] > 20)
                	rectangle(maskImg, Point(y, x), Point(y + (maskImg.cols / colDivisor) - 1, x + (maskImg.rows / rowDivisor) - 1), CV_RGB(255, 0, 0), 1); // visualization

                // imshow("Image", maskImg); // visualization
                // waitKey(0); // visualization
                k++;
            }
        }
        imshow("Segmented", maskImg); // visualization
        waitKey(0); // visualization
    }else if(img.cols % colDivisor != 0)
    {
        cerr << "Error: Please use another divisor for the column split." << endl;
        exit(1);
    }else if(img.rows % rowDivisor != 0)
    {
        cerr << "Error: Please use another divisor for the row split." << endl;
        exit(1);
    }
return EXIT_SUCCESS;
}

int segment(Mat &block, float thres) {
	Scalar mean, stddev;
	meanStdDev(block, mean, stddev);
	if (stddev.val[0] < thres) {
		block.setTo(Scalar(0,0,0));
		return 0;
	}
	return 1;
}

int makemask(Mat &block, float thres) {
	Scalar mean, stddev;
	meanStdDev(block, mean, stddev);
	if (stddev.val[0] < thres) {
		block.setTo(Scalar(0,0,0));
	} else {
		block.setTo(Scalar(255,255,255));
	}
}

int calculate_orientation(Mat &block) {
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

  	// cout << abs_grad_y << endl;
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
  	cout << angle << endl;

	int length = 10;
	Point P1(0,0);
	if (angle < 0) {
		angle *= -1;
	} else {
		angle *= -1;
		P1.x = 0;
		P1.y = 16;
	}
	Point P2;

	P2.x =  (int)round(P1.x + length * cos(angle * CV_PI / 180.0));
	P2.y =  (int)round(P1.y + length * sin(angle * CV_PI / 180.0));

	cout << P1 << " " << P2 << endl;


		line(block, P1, P2, CV_RGB(0, 255, 0));
		imshow("Block", block);
		// waitKey(0);
		block.setTo(Scalar(0,0,0));
		line(block, P1, P2, CV_RGB(0, 255, 0));
	
}

int main(int argc, char** argv) {
	//Read image
	String img = argv[1];
	Mat image = imread(img);

	//Check for failure
	if (image.empty()) {
		cout << "Could not open or find the image" << endl;
		return -1;
	}

	// normalize(image, -1, 1,)
	Mat normim;

	String windowName = "Fingerprint"; //Name of window
	namedWindow(windowName); //Create a window
	imshow(windowName, image); //Show image inside window
	// waitKey(0); //Wait for keystroke in window
	// destroyWindow(windowName); //Destroy

	// vector<Mat> BGR;
	// split(image, BGR);
	// cout << int(BGR[1].at<uchar>(0,0)) << endl;
	// cout << image.at<Vec3b>(0,0) << endl;

	cv::Scalar mean, stddev;
    cv::meanStdDev(image, mean, stddev);
    cout << "Mean " << mean.val[0] << endl;
    cout << "Stddev " << stddev.val[0] << endl;

    subtract(image, 1.0f*mean.val[0], normim);
    // normim = normim*1.0f/stddev.val[0];
    imshow("Normim", normim);
    waitKey(0);

 
    // for (int i=0 ; i<normim.rows ; i++) {
    // 	for (int j=0 ; j<normim.cols ; j++) {
    // 		cout << normim.at<double>(i,j) << " ";
    // 	}
    // 	cout << endl;
    // }

    meanStdDev(normim, mean, stddev);
    cout << "Mean " << mean.val[0] << endl;
    cout << "Stddev " << stddev.val[0] << endl;

	// waitKey(0);

	std::vector<cv::Mat> blocks;

	subdivide(image, 16, 16, blocks);

	

	for (int i=0 ; i<256 ; i++) {
		// imshow("block", blocks[i]);
		// waitKey(0);
		mask[i] = segment(blocks[i], 20);
		// calculate_orientation(blocks[i]);
		// if (i == 128) {
		// 	// cout << blocks[i] << endl;
		// 	Point x(0,0), y(10,10);
		// 	// rectangle(blocks[i], x,y,CV_RGB(0,255,0));
		// 	imshow("Block", blocks[i]);
		// }

	}

	Mat segmented = image.clone(); 	

	// makemask(mask, 20);

	// imshow("Mask", mask);

	for (int i=0; i<256; i++) {
		if (mask[i]) calculate_orientation(blocks[i]);
	}

	imshow("Image", segmented);
	imshow("Orientation", image);
	waitKey(0);

	int angle =-60;
	int length = 8;
	Point P1(50,50);
	Point P2;

	P2.x =  (int)round(P1.x + length * cos(angle * CV_PI / 180.0));
	P2.y =  (int)round(P1.y + length * sin(angle * CV_PI / 180.0));

	line(image, P1, P2, CV_RGB(0, 0, 255));
	line(image, P1, Point(60,50), CV_RGB(255, 0, 0));
	// imshow("Lala", image);
	cout << P1 << " " << P2 << endl;
	waitKey(0);
	return 0;
}

//g++ -ggdb test.cpp -o test `pkg-config --cflags --libs opencv`

/*
int main() {
	Mat image(Size(16,16), CV_8U);

	for(int i=0; i<image.rows ; i++) {
		for(int j=0; j<image.cols ; j++) {
			image.at<uchar>(i,j) = i*16 + j;
			// cout << i*16 + j << endl;
		}
	}

	cout << image << endl;
	imshow("image", image);
	waitKey(0);

	cv::Scalar mean, stddev;
    cv::meanStdDev(image, mean, stddev);
    cout << "Mean " << mean.val[0] << endl;
    cout << "Stddev " << stddev.val[0] << endl;

	return 0;
}*/