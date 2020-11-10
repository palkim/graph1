#include <iostream>
#include <math.h> //I add this
#include "parser.h"
#include "ppm.h"

typedef unsigned char RGB[3];

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

parser::Vec3f cross(parser::Vec3f a, parser::Vec3f b) {
    
    parser::Vec3f tmp;

    tmp.x = a.y*b.z-b.y*a.z;
    tmp.y = b.x*a.z-a.x*b.z;
    tmp.z = a.x-b.y-b.x*a.y;

    return tmp;

}


int main(int argc, char* argv[])
{
    // Sample usage for reading an XML scene file
    parser::Scene scene;

    scene.loadFromXml(argv[1]);

    //general data
    parser::Vec3i background_color = scene.background_color;
    float shadow_ray_epsilon = scene.shadow_ray_epsilon;
    int max_recursion_depth = scene.max_recursion_depth;
    std::vector<parser::Camera> cameras = scene.cameras;
    parser::Vec3f ambient_light = scene.ambient_light;
    std::vector<parser::PointLight> point_lights = scene.point_lights;
    std::vector<parser::Material> materials = scene.materials;
    std::vector<parser::Vec3f> vertex_data = scene.vertex_data;
    std::vector<parser::Mesh> meshes = scene.meshes;
    std::vector<parser::Triangle> triangles = scene.triangles;
    std::vector<parser::Sphere> spheres = scene.spheres;

    //camera data
    parser::Vec3f position;
    parser::Vec3f gaze;
    parser::Vec3f up;
    parser::Vec3f parallel;
    parser::Vec4f near_plane;
    float near_distance;
    int image_width;
    int image_height;
    std::string image_name;
    double pixel_width;
    double pixel_height;
    double half_pixel_width;
    double half_pixel_height;


    for (int c = 0; c < cameras.size(); c++) {
        
        position = cameras[c].position;
        gaze = cameras[c].gaze;
        up = cameras[c].up;
        near_plane = cameras[c].near_plane; //left(x),right(y),bottom(z),top(w)
        near_distance = cameras[c].near_distance;
        image_width = cameras[c].image_width;
        image_height = cameras[c].image_height;
        image_name = cameras[c].image_name;

        parallel = cross(gaze, up);
        parallel = normalize(parallel);

        up = cross(parallel, gaze);
        up = normalize(up);

        gaze = normalize(gaze);
        
        pixel_width = (near_plane.x - near_plane.y)/(double) image_width;
        half_pixel_width = pixel_width*0.5;

        pixel_height = (near_plane.w - near_plane.z)/(double) image_height;
        half_pixel_height = pixel_height*0.5;

        char* image_name_ptr = &image_name[0];
        

        unsigned char* image = new unsigned char [image_width * image_height * 3];

        int i = 0;
        for (int y = 0; y < image_height; ++y)
        {
            for (int x = 0; x < image_width; ++x)
            {
                image[i++] = 155;
                image[i++] = 155;
                image[i++] = 155;
            }
        }

        write_ppm(image_name_ptr, image, image_width, image_height);
    }

    // The code below creates a test pattern and writes
    // it to a PPM file to demonstrate the usage of the
    // ppm_write function.
    //
    // Normally, you would be running your ray tracing
    // code here to produce the desired image.

    /*const RGB BAR_COLOR[8] =
    {
        { 255, 255, 255 },  // 100% White
        { 255, 255,   0 },  // Yellow
        {   0, 255, 255 },  // Cyan
        {   0, 255,   0 },  // Green
        { 255,   0, 255 },  // Magenta
        { 255,   0,   0 },  // Red
        {   0,   0, 255 },  // Blue
        {   0,   0,   0 },  // Black
    };

    //int width = 640, height = 480;
    //int columnWidth = width / 8;

    unsigned char* image = new unsigned char [width * height * 3];

    int i = 0;
    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; ++x)
        {
            int colIdx = x / columnWidth;
            image[i++] = BAR_COLOR[colIdx][0];
            image[i++] = BAR_COLOR[colIdx][1];
            image[i++] = BAR_COLOR[colIdx][2];
        }
    }

    write_ppm("test.ppm", image, width, height);
    */
}
