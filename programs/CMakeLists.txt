include_directories(.)

# Create programs.
add_executable(pm_pm4gnn pm4gnn.c)
add_executable(pm_dglpart dglpart.c)

# Link with the required libraries 
foreach(prog pm_pm4gnn pm_dglpart)
  target_link_libraries(${prog} pm4gnn metis GKlib m)
endforeach(prog)

install(TARGETS pm_pm4gnn pm_dglpart 
  RUNTIME DESTINATION bin)
