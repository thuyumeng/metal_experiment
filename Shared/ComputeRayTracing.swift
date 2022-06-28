//
//  ComputeRayTracing.swift
//  raytracing_tutorial_swift
//
//  Created by 于蒙 on 2022/6/21.
//

import Foundation
import Metal
import MetalKit


struct Vec3 {
    // 生成uniform分布的单位向量,用rejection method采样生成
    static func random_in_unit_sphere() ->Vec3{
        while(true){
            let v = Vec3(x: Float.random(in: -1.0..<1.0),
                         y: Float.random(in: -1.0..<1.0),
                         z: Float.random(in: -1.0..<1.0))
            if (v.length_squared() >= 1) {
                continue
            }
            return v
        }
    }
    
    var x: Float
    var y: Float
    var z: Float
    init(x: Float, y: Float, z:Float) {
        self.x = x
        self.y = y
        self.z = z
    }
    
    func length_squared() -> Float {
        return x*x + y*y + z*z
    }
}
 
func * (left: Float, right: Vec3) -> Vec3 {
    return Vec3(x: Float(left * right.x), y: left * right.y, z: left * right.z)
}
 
func + (left: Vec3, right: Vec3) -> Vec3 {
    return Vec3(x: left.x + right.x, y: left.y + right.y, z: left.z + right.z)
}
 
func - (left: Vec3, right: Vec3) -> Vec3 {
    return Vec3(x: left.x - right.x, y: left.y - right.y, z: left.z - right.z)
}
 
func dot (_ left: Vec3, _ right: Vec3) -> Float {
    return Float(left.x * right.x + left.y * right.y + left.z * right.z)
}
 
func unit_vector(_ v: Vec3) -> Vec3 {
    let length : Float = Float(sqrt(dot(v, v)))
    return Vec3(x: v.x/length, y: v.y/length, z: v.z/length)
}

struct Material {
    var material_type: UInt32
    var material_color: Vec3
    var fuzz: Float
    var ir: Float
}

struct metal_sphere {
    var center: Vec3
    var radius: Float
    var mtl: Material
}

struct metal_ray {
    var origin: Vec3
    var direction: Vec3
}

struct metal_camera {
    var lookfrom: Vec3
    var lookat: Vec3
    var vup: Vec3
    
    var vfov: Float
    var aspect_ratio: Float
}


func ComputeRayTracer(_ coordinator: Coordinator, _ view: MTKView){
    guard let drawable = view.currentDrawable else {
        return
    }
    
    let drawable_size = view.drawableSize
    let win_width = drawable_size.width
    let win_height = drawable_size.height
    
    do{
        let device = coordinator.metalDevice
        let metal_lib = device?.makeDefaultLibrary()
        let metal_func = metal_lib?.makeFunction(name: "compute_ray")
        let func_pso = try device?.makeComputePipelineState(function: metal_func!)
        
        // 设置command queue，command_buffer，compute_encoder
        let cmd_queue = coordinator.metalCommandQueue
        let cmd_buff = cmd_queue?.makeCommandBuffer()
        let compute_encoder = cmd_buff?.makeComputeCommandEncoder()
        // 设置compute encoder的pipeline state和输入buffer，或者texture
        // 即将compute shader和输入的数据绑定（笔者的理解)
        compute_encoder?.setComputePipelineState(func_pso!)
        
        let texture = drawable.texture
        compute_encoder?.setTexture(texture, index: 0)
        // 设置materialType
        let Diffuse: UInt32 = 0
        let Metal: UInt32 = 1
        let Dielectric: UInt32 = 2
        
        let mtl_ground = Material(material_type:Diffuse,
                                  material_color:Vec3(x: 0.8,
                                                      y: 0.3,
                                                      z: 0.3),
                                  fuzz:1.0,
                                  ir:1.0)
        let mtl_center = Material(material_type: Diffuse,
                                  material_color: Vec3(x: 0.8,
                                                       y: 0.8,
                                                       z: 0.8),
                                  fuzz:0.7,
                                  ir:1.0)
        let mtl_left = Material(material_type: Dielectric,
                                material_color: Vec3(x: 0.72,
                                                     y: 1.02,
                                                     z: 0.52),
                                fuzz:0.3,
                                ir:1.5)
        let mtl_right = Material(material_type: Metal,
                                 material_color: Vec3(x: 1.05,
                                                      y: 0.95,
                                                      z: 0.15),
                                 fuzz:0.0,
                                 ir:1.0)
        // 设置hitlist
        var spheres = [metal_sphere]()
        spheres.append(
            metal_sphere(center:Vec3(x:0,y:-100.5,z:-1.0),
                         radius:100.0,
                         mtl:mtl_ground)
        )
        spheres.append(
            metal_sphere(center:Vec3(x:0,y:0,z:-1.0),
                         radius:0.5,
                         mtl:mtl_center)
        )
        spheres.append(
            metal_sphere(center:Vec3(x:-1.0,y:0,z:-1.0),
                         radius:0.5,
                         mtl:mtl_left)
        )
        spheres.append(
            metal_sphere(center:Vec3(x:-1.0,y:0,z:-1.0),
                         radius:-0.4,
                         mtl:mtl_left)
        )
        spheres.append(
            metal_sphere(center:Vec3(x:1.0,y:0,z:-1.0),
                         radius:0.5,
                         mtl:mtl_right)
        )
        
        let array_size = MemoryLayout<metal_sphere>.size * spheres.count
        let sphere_buffer = device?.makeBuffer(
            bytes: spheres,
            length: array_size,
            options: MTLResourceOptions.storageModeShared)
        compute_encoder?.setBuffer(sphere_buffer, offset: 0, index: 0)
        
        var sphere_cnts = [Int]()
        sphere_cnts.append(spheres.count)
        
        let cnts_buffer = device?.makeBuffer(
            bytes: sphere_cnts,
            length: MemoryLayout<Int>.size,
            options: MTLResourceOptions.storageModeShared)
        compute_encoder?.setBuffer(cnts_buffer, offset:0, index:1)
        
        // Camera参数
        var camera_data = [metal_camera]()
        camera_data.append(
            metal_camera(lookfrom: Vec3(x:-2,
                                        y:2,
                                        z:1),
                         lookat: Vec3(x: 0,
                                      y: 0,
                                      z: -1),
                         vup: Vec3(x:0,
                                   y:1,
                                   z:0),
                         vfov: 40.0,
                         aspect_ratio: 16.0 / 9.0)
        )
        let buf_size = MemoryLayout<metal_camera>.size * camera_data.count
        
        let cam_buffer = device?.makeBuffer(
            bytes: camera_data,
            length: buf_size,
            options: MTLResourceOptions.storageModeShared)
        compute_encoder?.setBuffer(cam_buffer,
                                   offset: 0,
                                   index: 2)
        
        // 设置thread组织形式
        let grid_size = MTLSizeMake(Int(win_width), Int(win_height), 1)
        // metal 这个thread_group_size 应该怎样设置？
        // 参照 https://developer.apple.com/documentation/metal/calculating_threadgroup_and_grid_sizes
        // 这片文章：https://developer.apple.com/documentation/metal/creating_threads_and_threadgroups
        /*  介绍了 grid, threadgroup, simd_group三个层次的thread组织形式，对上面的文章是必不可少的说明补充*/
        let w = func_pso?.threadExecutionWidth
        let h = (func_pso?.maxTotalThreadsPerThreadgroup)! / w!
        let thread_group_size = MTLSizeMake(w!, h, 1)
        compute_encoder?.dispatchThreads(grid_size, threadsPerThreadgroup: thread_group_size)
        //  compute pass encoding结束
        compute_encoder?.endEncoding()
        
        cmd_buff?.present(drawable)
        // 执行command buffer
        cmd_buff?.commit()
    }
    catch {
        print("create function pipeline state failed")
    }
}
