add_library(enso-cap)

target_sources(enso-cap
  PRIVATE
    "src/enso/nonrigorous.cpp"
    "src/enso/rigorous.cpp"
    "src/enso/utility.cpp"
  PUBLIC
    FILE_SET HEADERS
    BASE_DIRS "include"
    FILES
      "include/enso-cap/enso/equation.h"
      "include/enso-cap/enso/nonrigorous.h"
      "include/enso-cap/enso/rigorous.h"
      "include/enso-cap/enso/utility.h"
      "include/enso-cap/enso-cap.h"
)
target_link_libraries(enso-cap PUBLIC capdDDEs::capdDDEs)
