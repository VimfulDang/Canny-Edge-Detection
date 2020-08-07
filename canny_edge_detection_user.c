#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <time.h>
//#include <intelfpgaup/video.h>
#define PI 3.141592654


/*
	1. Grayscale
	2. Gaussian Smoothing
	3. Sobel Operator (gradient)
	4. Non-maximum suppression
	5. Hysteresis

*/
typedef unsigned char byte;

char command[64];
int width, height;
int screen_x, screen_y, char_x, char_y;

struct pixel {
	byte b;
	byte g;
	byte r;
};


struct pixel zero_pixel = {
	.r = 0,
	.b = 0,
	.g = 0,
};

int read_bmp(char * filename, byte ** header, struct pixel ** data) {
	struct pixel *data_tmp;
	byte * header_tmp;
	int size;
	FILE * file = fopen (filename, "rb");

	if (!file) return -1;

	header_tmp = malloc (54 * sizeof(byte)); 

	fread (header_tmp, sizeof(byte), 54, file);

	//get height and width of image from the header
	width = *(int *) (header_tmp + 18); //18 byte offset
	height = *(int *) (header_tmp + 22); //22 byte offset

	// read in the image
	size = width * height;
	data_tmp = malloc (size * sizeof(struct pixel));
	fread(data_tmp, sizeof(struct pixel), size, file);

	*header = header_tmp;
	*data = data_tmp;
	return 0;
}

void write_bmp(char * filename, byte *header, byte * data) {
	FILE * file = fopen(filename, "wb");
	int y, x, size;

	//allow image[y] for a row  and image[y][x] for item in row
	struct pixel (*image)[width] = malloc(sizeof(struct pixel [height][width]));
	byte (*gray_image)[width] = (byte (*)[width]) data;
	fwrite(header, sizeof(byte), 54, file);

	//for grayscale images, copy blue and green value to be as red
	for (y = 0; y < height; y++) {
		for (x = 0; x < width; x++) {
			image[y][x].r = gray_image[y][x];
			image[y][x].b = gray_image[y][x];
			image[y][x].g = gray_image[y][x];
		}
	}

	size = width * height;
	//fwrite(source address, size of each item, no item, file descriptor)
	fwrite(image, sizeof(struct pixel), size, file);
	fclose(file);
	free(image);
}

// The input data is either the x- or y-derivative of the image, as calculated by Sobel. The
// output produced is a bmp file in which each pixel corresponds to the absolute value of the 
// derivative at that pixel. This bmp file allows us to visualize the derivative as a bmp image.
void write_signed_bmp(char *filename, byte *header, signed int *data) {
    FILE* file = fopen (filename, "wb");
    struct pixel bytes;
    int val;

    signed int (*image)[width] = (signed int (*)[width]) data; // allow image[][]
    // write the 54-byte header
    fwrite (header, sizeof(byte), 54, file); 
    int y, x;
    
    // convert the derivatives' values to pixels by copying each to an r, g, and b
    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            val = abs(image[y][x]);
            val = (val > 255) ? 255 : val;
            bytes.r = val;
            bytes.g = bytes.r;
            bytes.b = bytes.r;
            fwrite (&bytes, sizeof(struct pixel), 1, file); // write the data
        }
    }
    fclose (file);
}

//Convert image to grayscale by averaging rgb value and setting to same intensity.
void convert_to_grayscale(struct pixel * data, byte ** gray_ptr) {
	int x, y;
	struct pixel (*image)[width] = (struct pixel (*)[width]) data;
	byte (*gray_image)[width] = malloc (sizeof(byte [height][width]));

	for(y = 0; y < height; y++) {
		for (x = 0; x < width; x++) {
			gray_image[y][x] = (image[y][x].r + image[y][x].b + image[y][x].g) / 3;
		}
	}
	free(data);
	*gray_ptr = (byte *) gray_image;
}


// Render an image on the VGA Display
void draw_image(struct pixel * data, byte * gray_data, int fd_video, int grayscale) {
	int x, y, stride_x, stride_y, i, j, vga_x, vga_y;
	int r, g, b, color;
	struct pixel (*image)[width] = (struct pixel (*)[width]) data;
	byte (*image_gray)[width] = (byte (*)[width]) gray_data;

	//clear screen
	sprintf(command, "clear 0\n");
	write(fd_video, command, strlen(command));

	stride_x = (width > screen_x) ? width / screen_x : 1;
	stride_y = (height > screen_y) ? height / screen_y : 1;
	printf("stride_x:%d stride_y:%d\n", stride_x, stride_y);

	// scale proportionally (don't stretch image)

	stride_y = (stride_x > stride_y) ? stride_x : stride_y;
	stride_x = (stride_y > stride_x) ? stride_y : stride_x;

	for (y = 0; y < height; y += stride_y) {
		for (x = 0; x < width; x += stride_x) {
			if (grayscale) {
				r = 0;
				for (i = 0; i < stride_y; i++) {
					for(j = 0; j < stride_x; j++) {		
						r += image_gray[y][x];
					}
				}
				r = r / (stride_x * stride_y);
				color = (r >> 3) << 11 | (r >> 2) << 5 | (r >> 3);
			}
			else {
				//average colors for downsizing pixel
				r = 0, g = 0; b = 0;
				for (i = 0; i < stride_y; i++) {
					for(j = 0; j < stride_x; j++) {
						//first index is height's and second is within the array for width
						r += image[y + i][x + j].r;
						g += image[y + i][x + j].g;
						b += image[y + i][x + j].b;
					}
				}
				r /= (stride_x * stride_y);
				b /= (stride_x * stride_y);
				g /= (stride_x * stride_y);

				r >>= 3;
				g >>= 2;
				b >>= 3;
				color = r << 11 | g << 5 | b;
			}
			vga_x = x / stride_x;
			vga_y = y / stride_y;

            if (screen_x > width / stride_x)   // center if needed
                sprintf(command, "pixel %d %d %#x\n", vga_x + (screen_x-(width/stride_x))/2, (screen_y-1) - vga_y, color);
            else 
                if ((vga_x < screen_x) && (vga_y < screen_y))
                    sprintf(command, "pixel %d %d %#x\n", vga_x, (screen_y-1) - vga_y, color);       
        	write(fd_video, command, strlen(command));
		}
	}
}

int valid_pixel(int x, int y, int pad) {
	int valid = 1;
	if ((y < pad) | (y > height) | (x < pad) | (x > width)) 
		valid = 0;
	return valid;
}

void gaussian_blur_grayscale(byte **data) {

	int x, y, i, j;
	int pad = 2;
	unsigned long pixel_value;
	unsigned int filter[5][5] = {
        { 1,  4,  7,  4, 1 },
        { 4, 16, 26, 16, 4 },
        { 7, 26, 41, 26, 7 },
        { 4, 16, 26, 16, 4 },
        { 1,  4,  7,  4, 1 }
    };

    //declare an image[][] variable to access pixel data
    byte (*image)[width] = (byte (*)[width]) *data;
    // allocate memory for the 2-D covolution image
    byte (*convolution)[width] = malloc(sizeof(byte [height][width]));
    for(y = 0; y < height; y++) {
    	for (x = 0; x < width; x++) {
    		pixel_value = 0;
    		for(i = 0; i < 5; i++) {
    			for(j = 0; j < 5; j++) {
		    		if (valid_pixel(x+j-pad, y+i-pad, pad))
	    				pixel_value += ((unsigned int) image[y+i-pad][x+j-pad]) * filter[i][j];
    			}
    		}
    		convolution[y][x] = (byte) (pixel_value / 273);	
    	}
    }

    free (*data);
    *data = (byte *) convolution;
}

/* This function finds the image gradient.
 * Input:
 *      data is a 2-D array of image pixels
 * Outputs: 
 *      data is used to return the intensity of the image gradient
 *      conv_x is used to return the convolution with the sobel_x operator (Df/dx)
 *      conv_y is used to return the convolution with the sobel_y operator (Df/dy)
 */


void sobel_filter(byte **data, signed int **conv_x, signed int **conv_y) {
   
   	unsigned int x, y, i, j;
   	signed int pixel_dx, pixel_dy, pixel_gradient;

    // Definition of Sobel filter in horizontal direction (highlights vertical edges)
    int sobel_x[3][3] = {
        { 1,  0,  -1 },
        { 2,  0,  -2 },
        { 1,  0,  -1 }
    };
    // Definition of Sobel filter in vertical direction (highlights horizontal edges)
    int sobel_y[3][3] = {
        { -1, -2, -1 },     // y == 0 (bottom) row
        {  0,  0,  0 },     // y == 1 (middle) row
        {  1,  2,  1 }      // y == 2 (top) row
    };

    signed int (*G_x)[width] = malloc (sizeof(signed int [height][width])); 
    signed int (*G_y)[width] = malloc (sizeof(signed int [height][width]));     
    byte (*gradient)[width] = malloc (sizeof(byte [height][width])); 
    // declare an image[][] varable to access the pixel data
    byte (*image)[width] = (byte (*)[width]) *data;
    
    for (y = 0; y < height; y++) {
    	for(x = 0; x < width; x++) {
    		pixel_dy = 0;
    		pixel_dx = 0;
    		pixel_gradient = 0;
    		for(i = 0; i < 3; i++) {
    			for(j = 0; j < 3; j++) {
    				if (valid_pixel(x+j-1, y+i-1, 1)) {
	    				pixel_dy += ((signed int) image[y+i-1][x+j-1]) * sobel_y[i][j];
    					pixel_dx += ((signed int) image[y+i-1][x+j-1]) * sobel_x[i][j];
    				}
    			}
    		}
    		pixel_gradient = (abs(pixel_dy) + abs(pixel_dx)) / 2;
    		G_x[y][x] = pixel_dx;
    		G_y[y][x] = pixel_dy;
    		gradient[y][x] = (pixel_gradient > 255) ? 255 : (byte) pixel_gradient;
    	}
    }

    free (*data);   // the previous image data is no longer needed
    *data = (byte *) gradient;
    *conv_x = (signed int *) G_x;
    *conv_y = (signed int *) G_y;
}


static int get_word(char * command, int offset, char * strbuff, int len) {
	int k;
	for(k = 0; k < len; k++, offset++) {
		if (*(command + offset) != ' ' && *(command + offset) != '\0') {
			*(strbuff + k) = *(command + offset);
		}
		else {
			*(strbuff + k)  = '\0';
			k = 9;
		}
	}
	return offset;
}

byte max_compare(byte back, byte mid, byte forward) {
	byte val = 0;
	if ((mid > back) & (mid >= forward))
		val = mid;
	return val;
}

/* Divide circle into 8 sections corresponding to 4 different edge directions
 * Each edge has a tolerance of +/- 22.5 degrees
 */
int determine_edge(signed int dx, signed int dy) {
	#define degree_per_rad 180/PI
	#define tolerance_degree 22.5
	int i, j, angle_index;
	double angle_degree;
	unsigned int edge_matrix[4][4] = {
		{7, 8, 15,  0}, //Vertical
		{1, 2,  9, 10}, // '\' diagonal
		{3, 4, 11, 12}, // Horizontal
		{5, 6, 13, 14} // '/' diagonal
	};
	double angle_rad = atan((double) dy/(double) dx);
	//printf("angle_rad:%f\n", angle_rad);
	if (angle_rad < 0)
		angle_rad += 2*PI;
	angle_degree = angle_rad * degree_per_rad / tolerance_degree;
	angle_index = (int) round(angle_degree);
	angle_index %= 16;
	for(i = 0; i < 4; i++) {
		for(j = 0; j < 4; j++) {
			if (angle_index == edge_matrix[i][j])
			{
				return i;
			}
		}
	}
	printf("dy:%d, dx:%d, angel_rad:%f, angle_index:%d\n", dy, dx, angle_rad, angle_index);
	return -1;
}

/* Makes edges thinner */
void non_max_suppress(byte **data, signed int *sobel_x, signed int *sobel_y) {
	int x, y, edge;
	byte back, mid, forward;
    byte (*temp_data)[width] = malloc (sizeof(byte [height][width])); 
    byte (*image)[width] = (byte (*)[width]) *data;
    signed int (*G_x)[width] = (signed int (*)[width]) sobel_x; // enable G_x[][]
    signed int (*G_y)[width] = (signed int (*)[width]) sobel_y; // enable G_y[][]
    for (y = 1; y < height-1; y++) {
    	for (x = 1; x < width-1; x++) {
    			edge = determine_edge(G_x[y][x], G_y[y][x]);
    			mid = image[y][x];
    			switch (edge) {
    				//Vertical
    				case (0) : back = image[y][x-1];
    						   forward = image[y][x+1];
    						   break;
				    // '\' diagonal
    				case (1) : back = image[y-1][x+1];
    						   forward = image[y+1][x-1];
    						   break;
    				// Horizontal
    				case (2) : back = image[y-1][x];
    						   forward = image[y+1][x];
    						   break;
    				// '/' diagonal
    				case (3) : back = image[y-1][x-1];
    						   forward = image[y+1][x+1];
    						   break;
    				default : back = 0;
    						  forward = 0;
    			}
    			temp_data[y][x] = max_compare(back, mid, forward);
    	}
    }

    free(*data);
    *data = (byte *) temp_data;
}

// Only keep pixels that are next to at least one strong pixel.
void hysteresis_filter(byte ** data) {
	#define strong_pixel 50 // example value

	int x, y, i, j, valid;

	byte (*temp_data)[width] = malloc(sizeof(byte [height][width]));
	//pointer to array size width
	byte (*image)[width] = (byte (*)[width]) *data;

    for(y = 0; y < height; y++) {
    	for (x = 0; x < width; x++) {
    		valid = 0;
    		if (image[y][x] > strong_pixel) {
	    		for(i = 0; i < 2; i++) {
	    			for(j = 0; j < 2; j++) {
	    				if (valid_pixel(x+j, y+i, 1)) {
	    					if (image[y+i][x+j] > strong_pixel) {
	    						valid = 1;
	    						break;
	    					}
	    				}
	    			}
	    			if (valid) {
	    				break;
	    			}
	    		}
	    	}
	    	temp_data[y][x] = (valid) ? image[y][x] : 0;
	    }
	}

	free (*data);
	*data = (byte *) temp_data;
}

void get_resolution(int fd) {
	char buffer[256];
	char xStr[8], yStr[8];
	int offset = 0, i = 0, readNo;
	char readIn;

	while ((readNo = read(fd, &readIn, 1)) > 0) {
		buffer[i++] = readIn;
	}

	buffer[i] = '\0';
	printf("%s", buffer);
	while (buffer[offset++] != ':') {
		continue;
	}
	offset = get_word(buffer, offset, yStr, 8);
	while (buffer[offset++] != ':') {
		continue;
	}
	offset = get_word(buffer, offset, xStr, 8);

	screen_x = atoi(xStr);
	screen_y = atoi(yStr);
}

int main(int argc, char * argv[]) {

	struct pixel *image;
	byte * gray_image;
	byte * header;
	signed int *G_x, *G_y;
	int debug = 0, video = 0;
	time_t start, end;
	int fd_video;

	if (argc < 3) {
        printf("Usage: part1 [-d] [-v] <BMP filename>\n");
        printf("-d: produces debug output for each stage\n");
        printf("-v: draws the input and output images on a video-out display\n");
        return 0;
	}

	if ((fd_video = open("/dev/video", O_RDWR)) == -1) {
		printf("Error opening /dev/video: %s\n", strerror(errno));
		return -1;
	}

	debug = atoi(argv[2]);
	//get resolution
	get_resolution(fd_video);
	printf("x:%d y:%d\n", screen_x, screen_y);	
	//test clear
	sprintf(command, "clear 1\n");
	write(fd_video, command, strlen(command));
	printf("%s\n", command);
	read_bmp(argv[1], &header, &image);


	// Start measuring time
    start = clock ();


	//draw_image(image, gray_image, fd_video, 0);
	//getchar();
	convert_to_grayscale(image, &gray_image);
	if (debug) write_bmp ("stage0_grayscale.bmp", header, gray_image);
	//draw_image(image, gray_image, fd_video, 1);
	//getchar();
	gaussian_blur_grayscale(&gray_image);
	//draw_image(image, gray_image, fd_video, 1);
	if (debug) write_bmp ("stage1_gaussian.bmp", header, gray_image);
	//getchar();
	sobel_filter(&gray_image, &G_x, &G_y);
	//draw_image(image, gray_image, fd_video, 1);
	if (debug) {
        write_signed_bmp ("stage2_gradient_x.bmp", header, G_x);
        write_signed_bmp ("stage2_gradient_y.bmp", header, G_y);
        write_bmp ("stage2_gradient.bmp", header, gray_image);
    }
	//getchar();
	non_max_suppress(&gray_image, G_x, G_y);
	//draw_image(image, gray_image, fd_video, 1);
	if (debug) write_bmp ("stage3_nonmax_suppression.bmp", header, gray_image);
	//hysteresis_filter(&gray_image);

	end = clock();
    printf("TIME ELAPSED: %.0f ms\n", ((double) (end - start)) * 1000 / CLOCKS_PER_SEC);
	

	draw_image(image, gray_image, fd_video, 1);
	if (debug) write_bmp ("stage4_hysteresis.bmp", header, gray_image);
	free(gray_image);
	free(header);
	close(fd_video);
	return 0;
}
