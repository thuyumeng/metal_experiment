//
//  ray_compute.metal
//  raytracing_tutorial_swift
//
//  Created by 于蒙 on 2022/6/21.
//

#include <metal_stdlib>
#include "random/random_header.metal"
using namespace metal;

// 下面的函数有些用到了引用：https://stackoverflow.com/questions/54665905/how-to-define-functions-with-a-referenced-parameter-in-metal-shader-language-exc
//  以上这片文章很好的解释了如何应用reference和地址空间


#define SAMPLES_PER_PIXEL 100
#define MAX_DEPTH 50
#define PI 3.1415926
#define Z_CORRECTION  0.0001


//Vec3 相关函数
struct Vec3{
    float x, y, z;
    Vec3(){
        x = 0.0;
        y = 0.0;
        z = 0.0;
    }
    Vec3(float x, float y, float z){
        this->x = x;
        this->y = y;
        this->z = z;
    }
    
    float length_squared(){
        return x*x + y*y + z*z;
    }
};

Vec3 operator+(thread const Vec3& u, thread const Vec3& v){
    return Vec3(u.x + v.x,
                u.y + v.y,
                u.z + v.z);
}

Vec3 operator-(thread const Vec3& u, thread const Vec3& v){
    return Vec3(u.x - v.x,
                u.y - v.y,
                u.z - v.z);
}

Vec3 operator*(float t, thread const Vec3& v) {
    return Vec3(t * v.x,
                t * v.y,
                t * v.z);
}

Vec3 operator*(thread const Vec3& v, float t) {
    return t * v;
}

Vec3 operator*(thread const Vec3& v, thread const Vec3& t) {
    return Vec3(
                v.x * t.x,
                v.y * t.y,
                v.z * t.z);
}

// 公共函数
float dot(thread const Vec3& u, thread const Vec3& v)
{
    return u.x*v.x + u.y*v.y + u.z*v.z;
}

Vec3 cross(thread const Vec3& u, thread const Vec3& v)
{
    float3 u1 = float3(u.x, u.y, u.z);
    float3 v1 = float3(v.x, v.y, v.z);
    float3 c = cross(u1, v1);
    return Vec3(c.x, c.y, c.z);
}
// 归一化Vec3
Vec3 unit_vector(thread const Vec3& direction)
{
    float vec_len = length(float3(direction.x,
                      direction.y,
                      direction.z));
    float rcp_len = 1.0 / vec_len;
    Vec3 ret(direction.x * rcp_len,
             direction.y * rcp_len,
             direction.z * rcp_len);
    return ret;
}

// 反射
Vec3 reflect(thread const Vec3& u, thread const Vec3& n) {
    float3 in_ray = float3(u.x, u.y, u.z);
    float3 normal = float3(n.x, n.y, n.z);
    float3 out_ray = reflect(in_ray, normal);
    return Vec3(out_ray.x, out_ray.y, out_ray.z);
//    return u - 2.0 * dot(u,n) * n;
}

// 折射
Vec3 refract(thread const Vec3& u, thread const Vec3& n, float ir)
{
    float3 in_ray = float3(u.x, u.y, u.z);
    float3 normal = float3(n.x, n.y, n.z);
    float3 out_ray = refract(in_ray, normal, ir);
    return Vec3(out_ray.x, out_ray.y, out_ray.z);
}

// 产生unit sphere 中的vector
Vec3 random_in_unit_sphere(thread pcg32_random_t* rng)
{
    float phi = 2.0 * PI * randomF(rng);
    float cosTheta = 2.0 * randomF(rng) - 1.0;
    float u = randomF(rng);

    float theta = acos(cosTheta);
    float r = pow(u, 1.0 / 3.0);

    float x = r * sin(theta) * cos(phi);
    float y = r * sin(theta) * sin(phi);
    float z = r * cos(theta);

    return Vec3(x, y, z);
}


// Ray光线数据结构，就是射线啦。。
struct Ray {
    Vec3 origin;
    Vec3 direction;

    Ray(thread const Vec3& origin, thread const Vec3& direction)
    {
        this->origin = origin;
        this->direction = direction;
    }
    
    Vec3 point_at_parameter(float t) const{
        return origin + t*direction;
    }
};

// Material实现用参数material_type(ENUM)和接口函数来控制不同材质的光线反射，折射性质
enum MaterialType{
    Diffuse=0,
    Metal,
    Dielectric,
};
    
struct Material {
    uint material_type;
    Vec3 material_color;
    // Metal材质的fuzz参数
    float fuzz;
    // Dielectric材质的refractin参数
    float ir;
};
    
// 实现两个模块
// 1、hitrecord：记录每次ray intersect的信息
// 2、用于检测和光线相交的Sphere。
// 3、用于光线相交检测的hittable list
struct HitRecord{
    Vec3 p;
    Vec3 normal;
    float t;
    Material mtl;
    bool front_face;
    
    HitRecord(){
        p = Vec3();
        normal = Vec3();
        t = 0.0;
        front_face = false;
    }
    
    void set_face_normal(thread const Ray& ray, thread const Vec3& outward_normal){
        this->front_face = (dot(ray.direction, outward_normal) < 0);
        
        if (front_face)
        {
            normal = outward_normal;
        }
        else{
            normal = -1.0 * outward_normal;
        }
    }
};
    
bool material_scatter(thread const Ray& ray_in, thread const HitRecord& hit_rec, thread Vec3& attenuation, thread Ray& scattered, thread pcg32_random_t* rng){
    uint material_type = hit_rec.mtl.material_type;
    switch(material_type)
    {
        case Diffuse:
        {
            Vec3 scatter_direction = hit_rec.normal + random_in_unit_sphere(rng);
            scattered = Ray(hit_rec.p + Z_CORRECTION * hit_rec.normal, scatter_direction);
            attenuation = hit_rec.mtl.material_color * attenuation;
            return true;
        }
        case Metal:
        {
            Vec3 reflected = unit_vector(reflect(ray_in.direction, hit_rec.normal)) + hit_rec.mtl.fuzz * random_in_unit_sphere(rng);
            scattered = Ray(hit_rec.p + Z_CORRECTION * hit_rec.normal, reflected);
            attenuation = hit_rec.mtl.material_color * attenuation;
            return (dot(scattered.direction, hit_rec.normal) > 0);
        }
        case Dielectric:
        {
            float ir = hit_rec.mtl.ir;
            float refraction_ratio = hit_rec.front_face ? 1.0 / ir : ir;
            
            // 全反升判断
            Vec3 unit_direction = unit_vector(ray_in.direction);
            Vec3 unit_normal = unit_vector(hit_rec.normal);
            float cos_theta = min(dot(-1.0 * unit_direction, unit_normal), 1.0);
            float sin_theta = sqrt(1.0 - cos_theta*cos_theta);
            bool can_reflect = sin_theta * refraction_ratio >= 1.0;
            
            // Schlick Approximation判断
            float r0 = (1.0 - refraction_ratio) / (1.0 + refraction_ratio);
            r0 = r0*r0;
            float reflectance = r0 + (1.0 - r0) * pow((1.0 - cos_theta), 5.0);
            
            Vec3 scattered_dir;
            Vec3 start_pt;
            if (can_reflect || reflectance > randomF(rng)){
                scattered_dir = reflect(ray_in.direction, hit_rec.normal);
                start_pt = hit_rec.p + Z_CORRECTION * hit_rec.normal;
            }
            else
            {
                start_pt = hit_rec.p - Z_CORRECTION * hit_rec.normal;
                scattered_dir = refract(unit_direction,
                                        unit_normal,
                                        refraction_ratio);
            }
            scattered = Ray(start_pt, scattered_dir);
            attenuation = hit_rec.mtl.material_color * attenuation;
            return true;
        }
        default:
        {
            return false;
        }
    }
}

// hittable objects
struct Sphere {
    Vec3 center;
    float radius;
    Material mtl;
    
    Sphere(thread const Vec3& center, float radius, Material mtl){
        this->center = center;
        this->radius = radius;
        this->mtl = mtl;
    }
    
    Sphere(device const Vec3& center, float radius, Material mtl){
        this->center = center;
        this->radius = radius;
        this->mtl = mtl;
    }
    
    bool hit(thread const Ray& ray, float t_min, float t_max, thread HitRecord& hit_record) const
    {
        Vec3 oc = ray.origin - center;
        float a = dot(ray.direction, ray.direction);
        float half_b = dot(oc, ray.direction);
        float c = oc.length_squared() - radius*radius;
        
        float discriminant = half_b*half_b - a*c;
        if (discriminant < 0){
            return false;
        }
        
        float sqrtd = sqrt(discriminant);
        float root = (-half_b - sqrtd) / a;
        if (root < t_min || t_max < root){
            root = (-half_b + sqrtd) / a;
            if (root < t_min || t_max < root){
                return false;
            }
        }
        
        hit_record.t = root;
        hit_record.p = ray.point_at_parameter(root);
        Vec3 outward_normal = 1.0 / radius * (hit_record.p - center);
        hit_record.set_face_normal(ray, outward_normal);
        return true;
    }
};

struct HittableList{
    device const Sphere* spheres;
    int sphere_cnts;
    
    HittableList(device const Sphere* spheres, int sphere_cnts){
        this->spheres = spheres;
        this->sphere_cnts = sphere_cnts;
    }
    
    bool hit(thread const Ray& ray, float t_min, float t_max, thread HitRecord& hit_record) const{
        float hit_anything = false;
        float closet_so_far = t_max;
        HitRecord temp_rec = HitRecord();
        
        for (int i=0; i<sphere_cnts; i++){
            Sphere sphere = Sphere(spheres[i].center,
                                   spheres[i].radius,
                                   spheres[i].mtl);
            if (sphere.hit(ray, t_min, closet_so_far, temp_rec)){
                hit_anything = true;
                closet_so_far = temp_rec.t;
                hit_record = temp_rec;
                hit_record.mtl = sphere.mtl;
            }
        }
        return hit_anything;
    }
};

// 为了减少传输buffer的大小设置CameraData类
struct CameraData{
    Vec3 lookfrom;
    Vec3 lookat;
    Vec3 vup;
    
    float vfov; // vertical field of view in degrees
    float apect_ratio;
};

// 构建Camera
struct Camera{
    // 配置参数
    CameraData static_data;
    
    // 中间计算的参数
    float theta;
    float focal_length;
    float viewport_height;
    float viewport_width;
    
    // 形成ray的Vec参数
    Vec3 origin;
    Vec3 lower_left_corner;
    Vec3 horizontal;
    Vec3 vertical;
    
    Camera(device const CameraData& cam_data){
        static_data = cam_data;
        
        theta = static_data.vfov / 180.0 * PI;
        focal_length = 1.0;
        float h = tan(0.5*theta);
        viewport_height = 2.0 * h;
        viewport_width = 1.0 * static_data.apect_ratio * viewport_height;
        
        Vec3 w = unit_vector(static_data.lookfrom - static_data.lookat);
        Vec3 u = unit_vector(cross(static_data.vup, w));
        Vec3 v = cross(w, u);
        
        horizontal = viewport_width * u;
        // -1.0 处理metal坐标方向
        vertical = -1.0 * viewport_height * v;
        origin = static_data.lookfrom;
        lower_left_corner = origin - 0.5*horizontal - 0.5*vertical - w;
    }
    
    Ray get_ray(float u, float v) {
        return Ray(origin, lower_left_corner + u*horizontal + v*vertical - origin);
    }
};

// ray_color 是反射天光的漫反射，每次反射回衰减
Vec3 ray_color(thread const Ray& ray, thread const HittableList& world, thread pcg32_random_t* rng)
{
    Ray cur_ray = ray;
    Vec3 cur_attenuation = Vec3(1.0, 1.0, 1.0);
    for (int i = 0; i<MAX_DEPTH; i++)
    {
        HitRecord rec;
        if(world.hit(cur_ray, 0.0, INFINITY, rec)){
            Ray scattered = cur_ray;
            if(material_scatter(cur_ray, rec, cur_attenuation, scattered, rng))
            {
                cur_ray = scattered;
            }
            else{
                break;
            }
        }
        else{
            Vec3 unit_direction = unit_vector(cur_ray.direction);
            float t = 0.5*(unit_direction.y + 1.0);
            Vec3 c = (1.0 - t)*Vec3(1.0, 1.0, 1.0) + t*Vec3(0.5, 0.7, 1.0);
            return cur_attenuation * c;
            
        }
    }
    return Vec3(0.0, 0.0, 0.0);
}

kernel void
compute_ray(device const Sphere* spheres,
            device const int* sphere_cnts,
            device const CameraData* cam_data,
            texture2d<half, access::write> outTexture [[texture(0)]],
            uint2                          gid        [[thread_position_in_grid]])
{
    // 设置camera常数
    Camera cam = Camera(cam_data[0]);
    // 多采样
    pcg32_random_t rng;
    uint64_t gidx = gid.x;
    uint64_t gidy = gid.y;
    pcg32_srandom_r(&rng, uint64_t(gidx *gidy), uint64_t(gidy));
    
    Vec3 color = Vec3();
    for(int i=0; i<SAMPLES_PER_PIXEL; i++)
    {
        float u = (float(gid.x) + randomF(&rng)) / float(outTexture.get_width() - 1);
        float v = (float(gid.y) + randomF(&rng)) / float(outTexture.get_height() - 1);
        Ray r = cam.get_ray(u, v);
        // 设置hittable list
        int cnt = sphere_cnts[0];
        HittableList world = HittableList(spheres, cnt);
        color = color + ray_color(r, world, &rng);
    }
    float rcp = 1.0 / float(SAMPLES_PER_PIXEL);
    color = rcp * color;

    float3 final_color = sqrt(float3(color.x,
                                     color.y,
                                     color.z));

    outTexture.write(half4(final_color.x, final_color.y, final_color.z, 1.0), gid);
    
}
