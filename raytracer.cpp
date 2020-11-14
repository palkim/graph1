#include <iostream>
#include <math.h> //I add this
#include <limits> // I add this
#include "parser.h"
#include "ppm.h"

typedef unsigned char RGB[3];

// data representing ray(t)=a(origin/camera cam_position)+t b(direction/calculated with gaze and distance between camera and image plane)
struct Ray
{
    parser::Vec3f a;
    parser::Vec3f b;
};
//data representing color
struct Color
{
    int R;
    int G;
    int B;
};

//general data
parser::Vec3i background_color;
float shadow_ray_epsilon;
int max_recursion_depth;
std::vector<parser::Camera> cameras;
parser::Vec3f ambient_light;
std::vector<parser::PointLight> point_lights;
std::vector<parser::Material> materials;
std::vector<parser::Vec3f> vertex_data;
std::vector<parser::Mesh> meshes;
std::vector<parser::Triangle> triangles;
std::vector<parser::Sphere> spheres;
int numSpheres;
int numTriangles;
int numMeshes;

//camera data
parser::Vec3f cam_position;
parser::Vec3f gaze;
parser::Vec3f up;
parser::Vec3f parallel;
parser::Vec4f near_plane;//left(x),right(y),bottom(z),top(w)
float near_distance;
int image_width;
int image_height;
std::string image_name;
double pixel_width;
double pixel_height;
double half_pixel_width;
double half_pixel_height;

char* image_name_ptr;
unsigned char* image;

double length(parser::Vec3f v) {
    return sqrt(v.x*v.x+v.y*v.y+v.z*v.z);
}

parser::Vec3f normalize(parser::Vec3f v) {
    parser::Vec3f tmp;
    double d;

    d = length(v);
    tmp.x = v.x/d;
    tmp.y = v.y/d;
    tmp.z = v.z/d;

    return tmp;
}

parser::Vec3f cross_product (parser::Vec3f a, parser::Vec3f b) {
    parser::Vec3f result;
    result.x = a.y * b.z - a.z * b.y;
    result.y = a.z * b.x - a.x * b.z;
    result.z = a.x * b.y - a.y * b.x;

    return result;
}

float dot_product (parser::Vec3f a, parser::Vec3f b) {
    float product = 0;
    product += a.x * b.x;
    product += a.y * b.y;
    product += a.z * b.z;

    return product;
}

parser::Vec3f scalar_multiplication(parser::Vec3f a, float b) {
    parser::Vec3f result;
    result.x = a.x * b;
    result.y = a.y * b;
    result.z = a.z * b;

    return result;
}

parser::Vec3f element_wise_multiplication(parser::Vec3f a, parser::Vec3f b) {
    parser::Vec3f result;
    result.x = a.x * b.x;
    result.y = a.y * b.y;
    result.z = a.z * b.z;

    return result;
}

parser::Vec3f element_wise_addition(parser::Vec3f a, parser::Vec3f b) {
    parser::Vec3f result;
    result.x = a.x + b.x;
    result.y = a.y + b.y;
    result.z = a.z + b.z;

    return result;
}

parser::Vec3f element_wise_subtraction(parser::Vec3f a, parser::Vec3f b) {
    parser::Vec3f result;
    result.x = a.x - b.x;
    result.y = a.y - b.y;
    result.z = a.z - b.z;

    return result;
}

void print_vector(parser::Vec3f a) {
    printf("(%f, %f, %f)", a.x, a.y, a.z);
    printf("\n");
}

void initImage() {
    int i = 0;
    for (int y = 0; y < image_height; ++y)
    {
        for (int x = 0; x < image_width; ++x)
        {
            image[i++] = background_color.x;
            image[i++] = background_color.y;
            image[i++] = background_color.z;
        }
    }
}

Color clamp(Color color) {
    Color result;
    
    if (color.R < 0) {
        result.R = 0;
    } else if (color.R > 250) {
        result.R = 250;
    } else {
        result.R = color.R;
    }

    if (color.G < 0) {
        result.G = 0;
    } else if (color.G > 250) {
        result.G = 250;
    } else {
        result.G = color.G;
    }

    if (color.B < 0) {
        result.B = 0;
    } else if (color.B > 250) {
        result.B = 250;
    } else {
        result.B = color.B;
    }

    return result;
} 

void readXml(char *fname) {
    parser::Scene scene;
    scene.loadFromXml(fname);

    //general data
    background_color = scene.background_color;
    shadow_ray_epsilon = scene.shadow_ray_epsilon;
    max_recursion_depth = scene.max_recursion_depth;
    cameras = scene.cameras;
    ambient_light = scene.ambient_light;
    point_lights = scene.point_lights;
    materials = scene.materials;
    vertex_data = scene.vertex_data;
    meshes = scene.meshes;
    triangles = scene.triangles;
    spheres = scene.spheres;

    numSpheres = spheres.size();
    numTriangles = triangles.size();
    numMeshes = meshes.size();
}

void setPixelInImage(int row, int col, Color color) {
    image[row*image_width*3+col*3] = color.R;
    image[row*image_width*3+col*3+1] = color.G;
    image[row*image_width*3+col*3+2] = color.B;
}

void setCameraData(parser::Camera cam) {
    cam_position = cam.position;
    gaze = cam.gaze;
    up = cam.up;
    near_plane = cam.near_plane; //left(x),right(y),bottom(z),top(w)
    near_distance = cam.near_distance;
    image_width = cam.image_width;
    image_height = cam.image_height;
    image_name = cam.image_name;

    parallel = cross_product(gaze, up);
    parallel = normalize(parallel);

    up = cross_product(parallel, gaze);
    up = normalize(up);

    gaze = normalize(gaze);
    
    pixel_width = (near_plane.x - near_plane.y)/(double) image_width;
    half_pixel_width = pixel_width*0.5;

    pixel_height = (near_plane.w - near_plane.z)/(double) image_height;
    half_pixel_height = pixel_height*0.5;
}

int main(int argc, char* argv[])
{
    if (argc < 2) {
        printf("Usage: ./raytracer <xml file>");
        return 1;
    }

    readXml(argv[1]);

    for (int c = 0; c < cameras.size(); c++) {
        setCameraData(cameras[c]);
        image_name_ptr = &image_name[0];

        image = new unsigned char [image_width*image_height*3];

        initImage();

        for (int row = 0; row < image_height; row++) {
            for (int col=0; col < image_width; col++) {
                ;
            }
        }
       
        write_ppm(image_name_ptr, image, image_width, image_height);
    }

    
    return 0;
}
