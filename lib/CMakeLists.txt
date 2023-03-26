# Include directories for library code.
include_directories(.)

# Find sources.
file(GLOB pm4gnn_sources *.c)

# Create libpm4gnn
add_library(pm4gnn ${PM4GNN_LIBRARY_TYPE} ${pm4gnn_sources})

install(TARGETS pm4gnn
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION lib
  ARCHIVE DESTINATION lib)
