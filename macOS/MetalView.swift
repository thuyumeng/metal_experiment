//
//  MetalView.swift
//  raytracing_tutorial_swift (iOS)
//
//  Created by 于蒙 on 2022/6/20.
//

import Foundation
import AVFoundation
import MetalKit
import SwiftUI


// 如果在swiftUI中加入MtkView需要更底层的UIKit和NSKit
// 分别对应iOS和macOS的ui开发sdk。
struct MetalView:NSViewRepresentable{
    
    func makeCoordinator() -> Coordinator {
        Coordinator(self)
    }
    
    func makeNSView(context: NSViewRepresentableContext<MetalView>) -> MTKView {
        let mtkView = MTKView()
        mtkView.delegate = context.coordinator
        mtkView.preferredFramesPerSecond = 30
        mtkView.enableSetNeedsDisplay = true
        if let metalDevice = MTLCreateSystemDefaultDevice() {
            mtkView.device = metalDevice
        }
        mtkView.framebufferOnly = false
        mtkView.clearColor = MTLClearColor(red: 0, green: 0, blue: 255, alpha: 0)
        mtkView.drawableSize = mtkView.frame.size
        return mtkView
    }
    
    func updateNSView(_ nsView: MTKView, context: NSViewRepresentableContext<MetalView>) {
    }
    
    
}

// 创建和swiftUI沟通的Cooridnator，名字和实现能将这个“协调器”
// 看的比较清楚，是双方交换信息的接口

class Coordinator: NSObject, MTKViewDelegate {
    init(_ parent: MetalView)
    {
        self.parent = parent
        if let metalDeivce = MTLCreateSystemDefaultDevice(){
            self.metalDevice = metalDeivce
        }
        self.metalCommandQueue = metalDevice.makeCommandQueue()
        
        super.init()
    }
    
    func update(_ view: MTKView)
    {
        print("check update coordinator")
    }
    
    func mtkView(_ view: MTKView, drawableSizeWillChange size: CGSize) {
    }
    
    func draw(in view: MTKView) {
        ComputeRayTracer(self, view)
    }

    
    var parent: MetalView
    var metalDevice: MTLDevice!
    
    var metalCommandQueue: MTLCommandQueue!
}
