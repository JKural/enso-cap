add_library(enso-cap)

target_sources(enso-cap
  PRIVATE
    "src/enso/nonrigorous.cpp"
    "src/enso/rigorous.cpp"
    "src/enso/utility.cpp"
    "src/enso-limit/nonrigorous.cpp"
    "src/enso-limit/rigorous.cpp"
    "src/enso-limit/utility.cpp"
)
# target_include_directories(enso-cap PUBLIC
#   # $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/enso-cap>
#   # $<INSTALL_INTERFACE:include/enso-cap>
#   ${CMAKE_CURRENT_SOURCE_DIR}/include/enso-cap
# )
target_include_directories(enso-cap
  PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:include>"
)
target_link_libraries(enso-cap PUBLIC capdDDEs::capdDDEs)
