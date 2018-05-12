#include <opencv2/opencv.hpp>
using namespace cv;
using namespace std;

int main(int argc, char** argv) {
	//Read image
	Mat image = imread("sao.jpg");

	//Check for failure
	if (image.empty()) {
		cout << "Could not open or find the image" << endl;
		return -1;
	}

	String windowName = "Fingerprint"; //Name of window
	namedWindow(windowName); //Create a window
	imshow(windowName, image); //Show image inside window
	waitKey(0); //Wait for keystroke in window
	destroyWindow(windowName); //Destroy

	vector<Mat> BGR;
	split(image, BGR);
	cout << int(BGR[1].at<uchar>(0,0)) << endl;
	cout << image.at<Vec3b>(0,0) << endl;

	waitKey(0);
	return 0;
}

//g++ -ggdb opencvtest.cpp -o opencvtest `pkg-config --cflags --libs opencv`