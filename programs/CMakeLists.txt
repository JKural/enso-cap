find_package(Boost REQUIRED)

function(add_simple_program name)
  add_executable(${name})
  target_sources(${name} PRIVATE "${name}.cpp")
  target_link_libraries(${name} PRIVATE enso-cap Boost::boost)
endfunction()

add_simple_program(find-fixed-point-candidate)
add_simple_program(fixed-point-proof)
add_simple_program(plot)
add_simple_program(searching-for-fixed-point-switch)
