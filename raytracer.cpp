#include <iostream>
#include <math.h> //I add this
#include <limits> // I add this
#include "parser.h"
#include "ppm.h"

typedef unsigned char RGB[3];


struct Ray
{
    parser::Vec3f origin ;
    parser::Vec3f direction ;
};

struct Color
{
    int R ;
    int G ;
    int B ;
};

struct ImagePlane
{
    float left ;
    float right ;
    float bottom ;
    float top ;
};

struct Object 
{
    int material_index ;
    parser::Vec3f normal_vector ;
    parser::Vec3f intersection_point ;
    bool isIntersects ;
    double intersects_at_t ;
};



//general data
parser::Vec3i background_color ;
float shadow_ray_epsilon ;
int max_recursion_depth ;
std::vector<parser::Camera> cameras ;
parser::Vec3f ambient_light ;
std::vector<parser::PointLight> point_lights ;
std::vector<parser::Material> materials ;
std::vector<parser::Vec3f> vertex_data ;
std::vector<parser::Mesh> meshes ;
std::vector<parser::Triangle> triangles ;
std::vector<parser::Sphere> spheres ;
int numSpheres ;
int numTriangles ;
int numMeshes ;

//camera data
parser::Vec3f e ;
parser::Vec3f w ; //w is opposite of gaze in slights, here it is gaze
parser::Vec3f v ;
parser::Vec3f u ;
ImagePlane image_plane ;
float camera_imagePlane_distance ;
int image_width ;
int image_height ;
std::string image_name ;
double pixel_width ;
double pixel_height ;
double half_pixel_width ;
double half_pixel_height ;

char* image_name_ptr ;
unsigned char* image ;

float length(parser::Vec3f a ) {
    return sqrt(a.x * a.x + a.y * a.y + a.z * a.z ) ;
}

parser::Vec3f normalize(parser::Vec3f a ) {
    parser::Vec3f result ;
    float d ;

    d = length(a) ;
    result.x = a.x/d ;
    result.y = a.y/d ;
    result.z = a.z/d ;

    return result ;
}

parser::Vec3f cross_product (parser::Vec3f a , parser::Vec3f b ) {
    parser::Vec3f result ;
    result.x = a.y * b.z - a.z * b.y ;
    result.y = a.z * b.x - a.x * b.z ;
    result.z = a.x * b.y - a.y * b.x ;

    return result ;
}

float dot_product (parser::Vec3f a , parser::Vec3f b ) {
    float product = 0 ;
    product += a.x * b.x ;
    product += a.y * b.y ;
    product += a.z * b.z ;

    return product ;
}

parser::Vec3f scalar_multiplication(parser::Vec3f a , float b ) {
    parser::Vec3f result ;
    result.x = a.x * b ;
    result.y = a.y * b ;
    result.z = a.z * b ;

    return result ;
}

parser::Vec3f element_wise_multiplication(parser::Vec3f a , parser::Vec3f b ) {
    parser::Vec3f result ;
    result.x = a.x * b.x ;
    result.y = a.y * b.y ;
    result.z = a.z * b.z ;

    return result ;
}

parser::Vec3f element_wise_addition(parser::Vec3f a , parser::Vec3f b ) {
    parser::Vec3f result ;
    result.x = a.x + b.x ;
    result.y = a.y + b.y ;
    result.z = a.z + b.z ;

    return result ;
}

parser::Vec3f element_wise_subtraction(parser::Vec3f a , parser::Vec3f b ) {
    parser::Vec3f result ;
    result.x = a.x - b.x ;
    result.y = a.y - b.y ;
    result.z = a.z - b.z ;

    return result ;
}

void print_vector(parser::Vec3f a ) {
    printf("(%f, %f, %f)", a.x, a.y, a.z ) ;
    printf("\n") ;
}

void setPixelInImage(int row , int col , Color color ) {
    image[row * image_width * 3 + col*3 ] = color.R ;
    image[row * image_width * 3 + col*3 + 1 ] = color.G ;
    image[row * image_width * 3 + col*3 + 2 ] = color.B ;
}

void initImage() {
    Color background = {background_color.x , background_color.y , background_color.z };
    for (int row = 0; row < image_height; row++ )
    {
        for (int col = 0; col < image_width; col++ )
        {
            setPixelInImage(row,col, background);
        }
    }
}

Color clamp(Color color ) {
    Color result ;
    
    if (color.R < 0 ) {
        result.R = 0 ;
    } else if (color.R > 250 ) {
        result.R = 250 ;
    } else {
        result.R = color.R ; 
    }

    if (color.G < 0 ) {
        result.G = 0 ;
    } else if (color.G > 250 ) {
        result.G = 250 ;
    } else {
        result.G = color.G ;
    }

    if (color.B < 0 ) {
        result.B = 0 ;
    } else if (color.B > 250 ) {
        result.B = 250 ;
    } else {
        result.B = color.B ;
    }

    return result ;
} 

void readXml(char *fname ) {
    parser::Scene scene ;
    scene.loadFromXml(fname ) ;

    background_color = scene.background_color ;
    
    shadow_ray_epsilon = scene.shadow_ray_epsilon ;
    max_recursion_depth = scene.max_recursion_depth ;
    
    cameras = scene.cameras ;
    
    ambient_light = scene.ambient_light ;
    point_lights = scene.point_lights ;
    
    materials = scene.materials ;
    vertex_data = scene.vertex_data ;
    
    meshes = scene.meshes ;
    triangles = scene.triangles ;
    spheres = scene.spheres ;

    numSpheres = spheres.size() ;
    numTriangles = triangles.size() ;
    numMeshes = meshes.size() ;
}

parser::Vec3f get_ray_point_at_t(Ray ray , float t ) {
    parser::Vec3f result;
    
    parser::Vec3f d_mult_t;
    d_mult_t = scalar_multiplication(ray.direction , t ) ;

    result = element_wise_addition(ray.origin , d_mult_t ) ;

    return result;
}

void setCameraData(parser::Camera cam ) {
    e = cam.position ;
    w = cam.gaze ;
    v = cam.up ;
    //near_plane = cam.near_plane ; left(x),right(y),bottom(z),top(w)
    image_plane.left = cam.near_plane.x ;
    image_plane.right = cam.near_plane.y ;
    image_plane.bottom = cam.near_plane.z ;
    image_plane.top = cam.near_plane.w ;
    camera_imagePlane_distance = cam.near_distance ;
    image_width = cam.image_width ;
    image_height = cam.image_height ;
    image_name = cam.image_name ;

    u = cross_product(w, v ) ;
    u = normalize(u ) ;

    v = cross_product(u , w ) ;
    v = normalize(v ) ;

    w = normalize(w ) ;
    
    pixel_width = (image_plane.left - image_plane.right ) / (float ) image_width ;
    half_pixel_width = pixel_width * 0.5 ;

    pixel_height = (image_plane.top - image_plane.bottom ) / (float ) image_height ;
    half_pixel_height = pixel_height * 0.5 ;
}

Ray get_ray_to_pixel (int row, int col ) {
    Ray result;
    
    parser::Vec3f m ;
    parser::Vec3f q ;
    parser::Vec3f s ;

    //construct m vector
    parser::Vec3f w_mult_distance ;
    w_mult_distance = scalar_multiplication(w, camera_imagePlane_distance ) ;
    m = element_wise_addition(e, w_mult_distance );

    //construct q vector
    parser::Vec3f lu ;
    parser::Vec3f tv ;
    parser::Vec3f lu_plus_tv;
    lu = scalar_multiplication(u, image_plane.left ) ;
    tv = scalar_multiplication(v, image_plane.top ) ;
    lu_plus_tv = element_wise_addition(lu, tv ) ;
    q = element_wise_addition(m, lu_plus_tv ) ;

    //construct s vector
    float su ;
    float sv ;
    parser::Vec3f u_mult_su ;
    parser::Vec3f v_mult_sv ;
    parser::Vec3f suu_subt_svv ;
    su = (col + 0.5 ) * (image_plane.right - image_plane.left ) / image_width ;
    sv = (row + 0.5 ) * (image_plane.top - image_plane.bottom ) / image_height ;
    u_mult_su = scalar_multiplication(u , su ) ;
    v_mult_sv = scalar_multiplication(v, sv ) ;
    suu_subt_svv = element_wise_subtraction(u_mult_su, v_mult_sv ) ;
    s = element_wise_addition(q, suu_subt_svv ) ;

    //construct result
    parser::Vec3f d ;
    d = element_wise_subtraction(s , e ) ;
    result = {e , d } ;

    return result;
}

Object intersect_ray_with_object (Ray ray ) {

    Object result ;
    result.isIntersects = false ;
    result.intersects_at_t = std::numeric_limits<double>::max() ;

    for (int sphere_index = 0 ; sphere_index < numSpheres ; sphere_index++ ) {
        parser::Sphere sphere = spheres[sphere_index] ;
        int sphere_m_id = sphere.material_id ;
        int sphere_c_id = sphere.center_vertex_id ;
        float sphere_r = sphere.radius ;

        parser::Material spere_m = materials[sphere_m_id] ;
        parser::Vec3f sphere_c = vertex_data[sphere_c_id] ;

        

    }
}

Color get_color_of_pixel (Ray ray , int recursion_depth ) {

    Object intersection_object = intersect_ray_with_object(ray) ;
    Color background = {background_color.x, background_color.y, background_color.z} ;
    if (!intersection_object.isIntersects) return background ;
    else {
        Color orange = {255, 165, 0} ;
        return orange ;
    }
}

int main(int argc , char* argv[] )
{
    if (argc < 2 ) {
        printf("Usage: ./raytracer <xml file>" );
        return 1 ;
    }

    readXml(argv[1]) ;

    for (int c = 0 ; c < cameras.size() ; c++ ) {
        setCameraData(cameras[c]) ;
        image_name_ptr = &image_name[0] ;

        image = new unsigned char [image_width*image_height*3 ] ;

        initImage() ;

        for (int row = 0 ; row < image_height ; row++ ) {
            for (int col=0 ; col < image_width ; col++ ) {
                
                Ray ray_to_pixel ;
                Color pixel_color ;

                ray_to_pixel = get_ray_to_pixel(row, col ) ;
                //below here recursion depth is zero since this is not a mirror ray
                pixel_color = get_color_of_pixel(ray_to_pixel, 0 ) ;

                setPixelInImage(row, col , pixel_color );
            }
        }
       
        write_ppm(image_name_ptr , image , image_width , image_height ) ;
    }

    return 0 ;

}
