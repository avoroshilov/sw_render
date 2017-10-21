# Simple software rasterizer
This project uses shared Core project (via git submodules) - in case of offline project assembly, please manually sync the Core repository and paste files into the `core` folder.

## Description
Experiments with recreating graphics pipeline, rendering is performed via software rasterization, without using 3D graphics APIs, like OpenGL, D3D or Vulkan.

Implements:
* Sutherland-Hodgman triangle clipping
* Back-/front-face culling
* Depth buffer
* Vertex attributes: coordinates, normal, texture coordinates, color
* Perspective corrected rasterization
* Transformations stack (modelview and projection)
* Spherical Environment Mapping

Several debug switches present, see `render/source/drawer.cpp` (e.g. disable perspective correction, render normals as color, etc.).

## Screenshots
Simple rotating textured cube:

<img src="examples/cube.jpg" alt="Textured cube" />

Simple rotating textured cube, without perspective correction:

<img src="examples/cube_noperspcorr.jpg" alt="Textured cube (no perspective correction)" />

Realtime deformed sphere mesh with normals:

<img src="examples/normals.jpg" alt="Normals rendered via line primitives" />

Deformed sphere with Spherical Environment Mapping:

<img src="examples/sem.jpg" alt="Spherical Environment Mapping" />


## License
Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International Public License
