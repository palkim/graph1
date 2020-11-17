#include <iostream>
#include <math.h> //I add this
#include <limits> // I add this
#include "parser.h"
#include "ppm.h"
#include <ctime>

typedef unsigned char RGB[3];
double distance_check;

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
    double left ;
    double right ;
    double bottom ;
    double top ;
};

struct Object 
{
    parser::Material material;
    parser::Vec3f normal_vector ;
    parser::Vec3f intersection_point ;
    bool isIntersects ;
    double intersects_at_t ;
};



//general data
parser::Vec3i background_color ;
double shadow_ray_epsilon ;
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
int numPointLights ;

//camera data
parser::Vec3f e ;
parser::Vec3f w ; //w is opposite of gaze in slights, here it is gaze
parser::Vec3f v ;
parser::Vec3f u ;
ImagePlane image_plane ;
double camera_imagePlane_distance ;
int image_width ;
int image_height ;
std::string image_name ;

char* image_name_ptr ;
unsigned char* image ;

double length(parser::Vec3f a ) {
    return sqrt((double) a.x * (double) a.x + (double) a.y * (double) a.y + (double) a.z * (double) a.z ) ;
}

parser::Vec3f normalize(parser::Vec3f a ) {
    parser::Vec3f result ;
    double d ;

    d = (double) length(a) ;
    result.x = (double) a.x/d ;
    result.y = (double) a.y/d ;
    result.z = (double) a.z/d ;

    return result ;
}

parser::Vec3f cross_product (parser::Vec3f a , parser::Vec3f b ) {
    parser::Vec3f result ;
    result.x = (double) a.y * (double) b.z - (double) a.z * (double) b.y ;
    result.y = (double) a.z * (double) b.x - (double) a.x * (double) b.z ;
    result.z = (double) a.x * (double) b.y - (double) a.y * (double) b.x ;

    return result ;
}

double dot_product (parser::Vec3f a , parser::Vec3f b ) {
    double product = (double) 0 ;
    product += (double) a.x * (double) b.x ;
    product += (double) a.y * (double) b.y ;
    product += (double) a.z * (double) b.z ;

    return product ;
}

parser::Vec3f scalar_multiplication(parser::Vec3f a , double b ) {
    parser::Vec3f result ;
    result.x = a.x * (double) b ;
    result.y = a.y * (double) b ;
    result.z = a.z * (double) b ;

    return result ;
}

parser::Vec3f element_wise_multiplication(parser::Vec3f a , parser::Vec3f b ) {
    parser::Vec3f result ;
    result.x = (double) a.x * (double) b.x ;
    result.y = (double) a.y * (double) b.y ;
    result.z = (double) a.z * (double) b.z ;

    return result ;
}

parser::Vec3f element_wise_addition(parser::Vec3f a , parser::Vec3f b ) {
    parser::Vec3f result ;
    result.x = (double) a.x + (double) b.x ;
    result.y = (double) a.y + (double) b.y ;
    result.z = (double) a.z + (double) b.z ;

    return result ;
}

parser::Vec3f element_wise_subtraction(parser::Vec3f a , parser::Vec3f b ) {
    parser::Vec3f result ;
    result.x = (double) a.x - (double) b.x ;
    result.y = (double) a.y - (double) b.y ;
    result.z = (double) a.z - (double) b.z ;

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
    numPointLights = point_lights.size() ;
}

parser::Vec3f get_ray_point_at_t(Ray ray , double t ) {
    parser::Vec3f result;
    
    parser::Vec3f d_mult_t;
    d_mult_t = scalar_multiplication(ray.direction , (double) t ) ;

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
    double su ;
    double sv ;
    parser::Vec3f u_mult_su ;
    parser::Vec3f v_mult_sv ;
    parser::Vec3f suu_subt_svv ;
    su = ((double) col + 0.5 ) * ((double) image_plane.right - (double) image_plane.left ) / (double) image_width ;
    sv = ((double) row + 0.5 ) * ((double) image_plane.top - (double) image_plane.bottom ) / (double) image_height ;
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

double determinant(parser::Vec3f a, parser::Vec3f b, parser::Vec3f c){
    return a.x*(b.y*c.z - b.z*c.y) - b.x*(a.y*c.z - c.y*a.z) + c.x*(a.y*b.z - b.y*a.z);
}


Object intersect_ray_with_object (Ray ray , int row, int col, bool isShadow) {

    Object result ;
    result.isIntersects = false ;
    
    double intersects_t_max = std::numeric_limits<double>::max() ;

    for (int sphere_index = 0 ; sphere_index < numSpheres ; sphere_index++ ) {
        parser::Sphere sphere = spheres[sphere_index] ;
        int sphere_m_id = sphere.material_id ;
        int sphere_c_id = sphere.center_vertex_id ;
        parser::Material sphere_m ;
        parser::Vec3f sphere_c ; 
        parser::Vec3f ray_origin = ray.origin;
        parser::Vec3f ray_direction = ray.direction ;
        double sphere_r ;
        sphere_m = materials[sphere_m_id-1] ;
        
        sphere_c = vertex_data[sphere_c_id-1] ;
        sphere_r = sphere.radius ;

        double a ;
        double b ;
        double c ;
        double determinant ;
        
        a = dot_product(ray_direction , ray_direction ) ;  
        parser::Vec3f o_subt_c = element_wise_subtraction(ray_origin, sphere_c ) ;
        b = (double) 2 * dot_product(ray_direction , o_subt_c ) ;
        c = dot_product(o_subt_c , o_subt_c ) - (double) pow(sphere_r,2) ;
        determinant = (double) pow(b , 2 ) - ((double) 4 * a * c ) ;

        if (determinant < (double) 0) {
            //printf("determinant is zero, no intersection\n");
            continue ;
        } else {
            
            double t1 ;
            double t2 ;
            double eq1 = dot_product(ray_direction , ray_direction ) ;
            double eq2 = dot_product(scalar_multiplication(ray_direction, -1) , o_subt_c ) ;
            double determinant_sqrt = sqrt(determinant ) ;    

            t1 = (eq2 + determinant_sqrt ) / eq1 ;
            t2 = (eq2 - determinant_sqrt ) / eq1 ;
            //NOTE: should i put an intersection epsilon value ?
            
            if (t2 < (double) 0 && t1 < (double) 0) {
                continue ;
            } else {
                //printf("determinant is bigger than or equal to zero and there is intersection\n");
                result.isIntersects = true ;
                if (t1 < (double) 0 && t2 < intersects_t_max) {
                    intersects_t_max = t2 ;
                    result.intersects_at_t = t2 ;
                } else if (t2 < (double) 0 && t1 < intersects_t_max) {
                    intersects_t_max = t1 ;
                    result.intersects_at_t = t1 ;
                } else if (t1 < t2 && t1 < intersects_t_max) {
                    intersects_t_max = t1 ;
                    result.intersects_at_t = t1 ;
                } else if (t2 <= t1 && t2 < intersects_t_max) {
                    intersects_t_max = t2 ;
                    result.intersects_at_t = t2 ;
                }
            }
            result.material = sphere_m ;
            result.intersection_point = get_ray_point_at_t(ray, result.intersects_at_t) ;
            parser::Vec3f normal = element_wise_subtraction(result.intersection_point , sphere_c ) ;
            result.normal_vector = normalize(normal) ;

        }
    }
    parser::Face face;
    parser::Material triangle_m;
    
    parser::Vec3f a ;
    parser::Vec3f b  ;
    parser::Vec3f c ;

    parser::Vec3f c_subt_b  ;
    parser::Vec3f a_subt_b  ;

    parser::Vec3f normal ;

    parser::Vec3f beta_0;
    parser::Vec3f beta_1;
    parser::Vec3f beta_2;

    parser::Vec3f gama_0 ;
    parser::Vec3f gama_1 ;
    parser::Vec3f gama_2 ;

    parser::Vec3f t_0 ; //a-b
    parser::Vec3f t_1; //a-c
    parser::Vec3f t_2 ; //a-o

    double detA ;

    double beta ;
    double gama ;
    double t ;
    for (int triangle_index = 0; triangle_index < numTriangles; triangle_index++) {
        face = triangles[triangle_index].indices ;
        triangle_m = materials[triangles[triangle_index].material_id - 1] ;
        
        a = vertex_data[face.v0_id-1] ;
        b = vertex_data[face.v1_id-1] ;
        c = vertex_data[face.v2_id-1] ;

        c_subt_b = element_wise_subtraction(c , b) ;
        a_subt_b = element_wise_subtraction(a , b) ;

        normal = normalize(cross_product(c_subt_b, a_subt_b) ) ;

        beta_0 = element_wise_subtraction(a, ray.origin);
        beta_1 = element_wise_subtraction(a, c);
        beta_2 = ray.direction;

        gama_0 = element_wise_subtraction(a, b);
        gama_1 = element_wise_subtraction(a, ray.origin);
        gama_2 = ray.direction;

        t_0 = gama_0; //a-b
        t_1 = beta_1; //a-c
        t_2 = gama_1; //a-o


        detA = determinant(element_wise_subtraction(a, b), element_wise_subtraction(a, c), ray.direction);

        beta = determinant(beta_0, beta_1, beta_2) / detA;
        gama = determinant(gama_0, gama_1, gama_2) / detA;
        t = determinant(t_0, t_1, t_2) / detA;
        
        if(t < distance_check) {
            //printf("t < distance_check, no intersection") ;
            continue;
        }


        if((beta > (double) 0 || abs(beta) - (double)0 < 0.0001) && (gama > (double) 0 || abs(gama) - (double)0 < 0.0001) && ( ( beta + gama < (double) 1 ) || ( (beta + gama) - (double) 1 < 0.0001 ) ) ){

            if (t < intersects_t_max) {
                //printf("intersection occured") ;
                intersects_t_max = t ;
                result.intersects_at_t = t ;
                result.isIntersects = true ;
                result.intersection_point = get_ray_point_at_t(ray, t) ;
                result.material = triangle_m ;
                result.normal_vector = normal ; 
            }
        }
        
    }
    
    for (int mesh_index = 0 ; mesh_index < numMeshes ; mesh_index++) {
        parser::Mesh mesh = meshes[mesh_index] ;
        int numFaces = (mesh.faces).size() ;
        for (int face_index = 0; face_index < numFaces; face_index++) {
            face = mesh.faces[face_index] ;
            triangle_m = materials[mesh.material_id - 1] ;
            
            a = vertex_data[face.v0_id-1] ;
            b = vertex_data[face.v1_id-1] ;
            c = vertex_data[face.v2_id-1] ;

            c_subt_b = element_wise_subtraction(c , b) ;
            a_subt_b = element_wise_subtraction(a , b) ;

            normal = normalize(cross_product(c_subt_b, a_subt_b) ) ;

            beta_0 = element_wise_subtraction(a, ray.origin);
            beta_1 = element_wise_subtraction(a, c);
            beta_2 = ray.direction;

            gama_0 = element_wise_subtraction(a, b);
            gama_1 = element_wise_subtraction(a, ray.origin);
            gama_2 = ray.direction;

            t_0 = gama_0; //a-b
            t_1 = beta_1; //a-c
            t_2 = gama_1; //a-o


            detA = determinant(element_wise_subtraction(a, b), element_wise_subtraction(a, c), ray.direction);

            beta = determinant(beta_0, beta_1, beta_2) / detA;
            gama = determinant(gama_0, gama_1, gama_2) / detA;
            t = determinant(t_0, t_1, t_2) / detA;  

            if (row == 610 && col == 750 ) {
                printf("t: %f\n", t) ;
                printf("distance_check: %f\n", distance_check) ;
                printf("distance_check_val: %f\n", abs(t-distance_check));
                printf("distance_check_is_true: %d\n", abs(t-distance_check) < 0.00001 || t < distance_check);
            }

            if(abs(t - distance_check) < 0.00001 || t < distance_check) {
                continue;
            }
            /*if (row == 610 && col == 750) {
                printf("beta: %f\n", beta) ;
                printf("is beta bigger than zero: %d\n", beta > (double) 0) ;
                printf("is beta equal to zero: %d\n", (abs(abs(beta) - (double) 0)) < 0.0001) ;
                printf("gama: %f\n", gama) ;
                printf("is gama bigger than zero: %d\n", gama > (double) 0) ;
                printf("if gama equal to zero: %d\n", (abs(abs(gama) - (double) 0)) < 0.0001) ;
                printf("beta+gama: %f\n", beta+gama);
                printf("is beta+gama less than 1: %d\n", beta + gama < (double) 1);
                printf("is beta+gama is equal to 1: %d\n", (abs(beta+gama - (double) 1)) < 0.01) ;
            }*/
            
            if (beta > (double) 0 || abs(beta - (double) 0) < 0.001 ) {
                if (abs(beta - (double) 0) < 0.0001) {
                    beta = (double) 0 ;
                }
                if (gama > (double) 0 || abs(gama - (double) 0 ) < 0.001) {
                    if (abs(gama - (double) 0 ) < 0.0001) {
                        gama = (double) 0 ;
                    }
                    if (beta + gama < (double) 1 || (beta + gama - (double) 1) < 0.001) {
                        if (t < intersects_t_max) {
                            intersects_t_max = t ;
                            result.intersects_at_t = t ;
                            result.isIntersects = true ;
                            result.intersection_point = get_ray_point_at_t(ray, t) ;
                            result.material = triangle_m ;
                            result.normal_vector = normal ; 
                        }
                    }
                }
            }
            /*if((beta > (double) 0 || abs(abs(beta) - (double)0) < 0.015) && (gama > (double) 0 || abs(abs(gama) - (double)0) < 0.015) && ( ( beta + gama < (double) 1 ) || ( abs(abs(beta + gama) - (double) 1) < 0.015 ) ) ){
                
                if (t < intersects_t_max) {
                    intersects_t_max = t ;
                    result.intersects_at_t = t ;
                    result.isIntersects = true ;
                    result.intersection_point = get_ray_point_at_t(ray, t) ;
                    result.material = triangle_m ;
                    result.normal_vector = normal ; 
                }
            }*/
        }
    }
    if (row == 610 && col == 750) {
        printf("is shadow: %d\n", isShadow) ;
        printf("isIntersects: %d\n", result.isIntersects) ;
    }

    return result ;
}

Ray get_ray_for_shadow(parser::Vec3f intersection_point, parser::PointLight light) {
    Ray result ;

    parser::Vec3f wi = normalize(element_wise_subtraction(light.position, intersection_point)) ;
    parser::Vec3f wi_mult_eps = scalar_multiplication(wi, shadow_ray_epsilon) ;

    result.direction = wi ;    
    
    result.origin = element_wise_addition(intersection_point, wi_mult_eps) ;

    return result ;

}

parser::Vec3f spec_shading (parser::Vec3f point, parser::Vec3f normal, parser::Material material, parser::PointLight light, Ray ray) {
    parser::Vec3f result ;

    parser::Vec3f wi = element_wise_subtraction(light.position, point) ;
    wi = normalize(wi) ;

    parser::Vec3f normal_vector = normalize(normal) ;

    parser::Vec3f wo = element_wise_subtraction(ray.origin, point) ;
    wo = normalize(wo) ;

    parser::Vec3f wi_plus_wo = element_wise_addition(wi, wo) ;

    double one_over_wi_plus_wo_len = (((double) 1) / length(wi_plus_wo)) ;

    parser::Vec3f h = scalar_multiplication(wi_plus_wo, one_over_wi_plus_wo_len) ;

    double cos_theta = std::max((double) 0, dot_product(normal_vector, h)) ;

    double one_over_r_square = ((double) 1 ) / pow(length(element_wise_subtraction(light.position, point)),2) ;

    parser::Vec3f irradiance = scalar_multiplication(light.intensity, one_over_r_square) ;
    parser::Vec3f ks = material.specular ;
    result = element_wise_multiplication(ks, irradiance) ;
    result = scalar_multiplication(result, cos_theta) ;

    return result;

}

parser::Vec3f diff_shading (parser::Vec3f point, parser::Vec3f normal, parser::Material material, parser::PointLight light, Ray ray) {
    parser::Vec3f result ;

    parser::Vec3f wi = element_wise_subtraction(light.position, point) ;
    wi = normalize(wi) ;

    parser::Vec3f normal_vector = normalize(normal) ;

    parser::Vec3f wo = element_wise_subtraction(ray.origin, point) ;
    wo = normalize(wo) ;

    double cos_theta = std::max((double) 0, dot_product(wi, normal_vector)) ;

    double one_over_r_square = ((double) 1 ) / pow(length(element_wise_subtraction(light.position, point)),2) ;

    parser::Vec3f irradiance = scalar_multiplication(light.intensity, one_over_r_square) ;
    parser::Vec3f kd = material.diffuse ;
    result = element_wise_multiplication(kd, irradiance) ;
    result = scalar_multiplication(result, cos_theta) ;

    return result;

}

Color get_color_of_pixel (Ray ray , int recursion_depth , int row, int col) {

    Object intersection_object = intersect_ray_with_object(ray, row, col, false) ;
    Color background = {background_color.x, background_color.y, background_color.z} ;
    if (recursion_depth == 0) {
        distance_check = camera_imagePlane_distance;
    } else {
        distance_check = (double) 0;
    }
    if (!intersection_object.isIntersects) {
        return background ;
    } 
    else {
        Color pixel_color ;
        parser::Vec3f color_calc_temp = {0, 0, 0} ;
        parser::Material intersection_obj_mat = intersection_object.material ;
        parser::Vec3f intersection_point = intersection_object.intersection_point ;
        parser::PointLight light ;
        Ray shadow ;
        double light_point_dist ;
        Object shadow_intersection_object ;
        //ambient shading
        color_calc_temp = element_wise_addition(color_calc_temp, element_wise_multiplication(intersection_obj_mat.ambient, ambient_light)) ;
        
        //diffuse and specular shading
        for (int light_index = 0; light_index < numPointLights ; light_index++ ) {
            light = point_lights[light_index] ;
            shadow = get_ray_for_shadow(intersection_point, light) ;
            light_point_dist = length(element_wise_subtraction(light.position, intersection_point)) ;
            distance_check = (double) 0 ;
            shadow_intersection_object = intersect_ray_with_object(shadow, row, col, true) ;
            //printf("light_point_dist: %f\n", light_point_dist);
            //printf("shadow_int_at_t: %f\n", shadow_intersection_object.intersects_at_t);
            //printf("isIntersects: %d", shadow_intersection_object.isIntersects) ;
            //printf("\n") ;
            
            if (shadow_intersection_object.isIntersects && light_point_dist > shadow_intersection_object.intersects_at_t ) {
                //printf("dewam");
                continue;
            } else {
                color_calc_temp = element_wise_addition(color_calc_temp, diff_shading(intersection_point, intersection_object.normal_vector, intersection_obj_mat, light, ray)) ;
                color_calc_temp = element_wise_addition(color_calc_temp, spec_shading(intersection_point, intersection_object.normal_vector, intersection_obj_mat, light, ray)) ;
                //print_vector(shading(intersection_point, intersection_object.normal_vector, intersection_obj_mat, light, ray));            
            }
        }
        pixel_color = {round(color_calc_temp.x), round(color_calc_temp.y), round(color_calc_temp.z)} ;
        pixel_color = clamp(pixel_color) ;
        return pixel_color;
    }
}

int main(int argc , char* argv[] )
{
    clock_t begin = clock() ;
    if (argc < 2 ) {
        printf("Usage: ./raytracer <xml file>" );
        return 1 ;
    }

    readXml(argv[1]) ;

    for (int c = 0 ; c < cameras.size() ; c++ ) {
        if ( c != 2) {
            continue;
        }
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
                pixel_color = get_color_of_pixel(ray_to_pixel, 0 , row, col ) ;

                setPixelInImage(row, col , pixel_color );
            }
        }
        //Color orange = {255, 160, 0} ;
        //for (int p = 0 ; p < 50; p++) {
            
          //  setPixelInImage(610, 750+p, orange) ;
        //}
        //setPixelInImage(750, 800, orange) ;
        write_ppm(image_name_ptr , image , image_width , image_height ) ;
    }
    clock_t end = clock() ;
    double elapsed_secs = double(end-begin) / CLOCKS_PER_SEC ;
    printf("elapsed: %f", elapsed_secs);

    return 0 ;

}