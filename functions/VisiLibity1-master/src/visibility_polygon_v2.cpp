/**
 * \file visibility_polygon.cpp
 * \brief A MATLAB MEX file for computing a visibility polygon.
 * \author Andrew D. Horchler, horchler @ gmail.com
 * \date 2012-11-20
 * \version 1.0
 *
 * Copyright (c) 2012, Andrew D. Horchler
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * * Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 * * Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */
#include "visilibity.hpp" // VisiLibity header
#include "mex.h"          // MATLAB MEX header

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  // Check for proper number of arguments
  if (nrhs != 4) {
    mexErrMsgTxt("Four inputs required.");
  }
  if (nlhs > 2) { // Allow for a second output
    mexErrMsgTxt("Too many output arguments.");
  }
  if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
      mxGetNumberOfElements(prhs[0]) != 2) {
    mexErrMsgTxt("Input observer must be a noncomplex array of size 2.");
  }
  if (!mxIsCell(prhs[1])) {
    mexErrMsgTxt("Input environment must be a cell array.");
  }

  // Get observer coordinates
  double *observer_coords = mxGetPr(prhs[0]);
  VisiLibity::Point observer(observer_coords[0], observer_coords[1]);

  // Get environment polygons
  const mxArray *env_cell_array_ptr = prhs[1];
  mwSize num_polygons = mxGetNumberOfElements(env_cell_array_ptr);
  std::vector<VisiLibity::Polygon> polygons;
  for (int i = 0; i < num_polygons; i++) {
    const mxArray *poly_ptr = mxGetCell(env_cell_array_ptr, i);
    double *poly_coords = mxGetPr(poly_ptr);
    mwSize num_verts = mxGetM(poly_ptr);
    std::vector<VisiLibity::Point> vertices;
    for (int j = 0; j < num_verts; j++) {
      vertices.push_back(
          VisiLibity::Point(poly_coords[j], poly_coords[j + num_verts]));
    }
    polygons.push_back(VisiLibity::Polygon(vertices));
  }
  VisiLibity::Environment *my_environment =
      new VisiLibity::Environment(polygons);

  // Get epsilon and snap_distance
  double epsilon = mxGetScalar(prhs[2]);
  double snap_distance = mxGetScalar(prhs[3]);

  // Snap observer to the environment
  observer.snap_to_boundary_of((*my_environment), snap_distance);
  observer.snap_to_vertices_of((*my_environment), snap_distance);

  // Compute Visibility_Polygon
  VisiLibity::Visibility_Polygon my_vis_poly(observer, (*my_environment),
                                             epsilon);

  // Get the growing vertices from the visibility polygon object
  VisiLibity::Polygon growing_verts = my_vis_poly.get_growing_vertices();

  // Create an mxArray for the first output (the visibility polygon)
  plhs[0] = mxCreateDoubleMatrix(my_vis_poly.n(), 2, mxREAL);
  double *out_vis_poly = mxGetPr(plhs[0]);

  // Populate the first output array
  for (int i = 0; i < my_vis_poly.n(); i++) {
    out_vis_poly[i] = my_vis_poly[i].x();
    out_vis_poly[my_vis_poly.n() + i] = my_vis_poly[i].y();
  }

  // Check if the user has requested the second output argument (growing_verts)
  if (nlhs > 1) {
    // Create an mxArray for the second output
    plhs[1] = mxCreateDoubleMatrix(growing_verts.n(), 2, mxREAL);
    double *out_growing_verts = mxGetPr(plhs[1]);
    
    // Populate the second output array
    for (int i = 0; i < growing_verts.n(); i++) {
        out_growing_verts[i] = growing_verts[i].x();
        out_growing_verts[growing_verts.n() + i] = growing_verts[i].y();
    }
  }

  delete my_environment;

  return;
}