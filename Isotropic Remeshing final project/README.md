# Note
## Need to solve the search problem
## Need to complete the last part of Isotropic Remeshing

There seems to be an issue with igl::copyleft::cgal::extract_feature 
I had to type cast adj_faces[0] and adj_faces[1] on line 86 and 87 to 
int(adj_faces[0]) and int(adj_faces[1]) in the extract_feature.cpp of the libigl library

## After building, add the data folder into the build folder to use
