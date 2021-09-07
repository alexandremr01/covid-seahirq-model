solution "COVID-ITA-TASK"
   configurations { "Release", "Debug", "Optimized"}
   links {"Libraries"}
   defines {"NEA=16", "MAX_DAYS=400", "NA=13"}
   targetdir "../bin"
   includedirs { "../src",
                 "../3rdparty/eigen",
               }
   --buildoptions { "-std=c++11" }

   output_path = path.getabsolute("../output")
   output_path = "OUTPUT_PATH=\\\""..output_path.."\\\""
   output_path = output_path:gsub(" ", "\\ ")

   input_path = path.getabsolute("../input")
   input_path = "INPUT_PATH=\\\""..input_path.."\\\""
   input_path = input_path:gsub(" ", "\\ ")
   
   defines { output_path, input_path }

   ------------------------------ LIBS -----------------------------------------
   project "Libraries"
      kind "StaticLib"
      language "C++"
      files { "../src/solvers.cpp",
              "../src/real_data_reader.cpp",
              "../src/scenario.cpp",
              "../src/models.cpp",
            }

      configuration "Release"
         flags { "OptimizeSize", "OptimizeSpeed" }

      configuration "Debug"
         defines { "DEBUG" }
         flags { "Symbols" }

      configuration "Optimized"
         flags { "OptimizeSize", "OptimizeSpeed", "EnableSSE", "EnableSSE2" }

   ------------------------------ BINS -----------------------------------------
   project "spatial_covid0d_estrat"
      kind "ConsoleApp"
      language "C++"
      files { "../src/spatial_covid0d_estrat.cpp" }
      flags {"ExtraWarnings"}

      configuration "Release"
         flags { "OptimizeSize", "OptimizeSpeed" }

      configuration "Debug"
         defines { "DEBUG" }
         flags { "Symbols" }

      configuration "Optimized"
         flags { "OptimizeSize", "OptimizeSpeed", "EnableSSE", "EnableSSE2" }
   ----------------------------------------------------------------------------
   project "csv_to_input"
      kind "ConsoleApp"
      language "C++"
      files { "../src/csv_to_input.cpp" }
      flags {"ExtraWarnings"}

      configuration "Release"
         flags { "OptimizeSize", "OptimizeSpeed" }

      configuration "Debug"
         defines { "DEBUG" }
         flags { "Symbols" }

      configuration "Optimized"
         flags { "OptimizeSize", "OptimizeSpeed", "EnableSSE", "EnableSSE2" }
   ----------------------------------------------------------------------------
   -- project "spatial_covid0d_estrat2"
   --    kind "ConsoleApp"
   --    language "C++"
   --    files { "../src/spatial_covid0d_estrat2.cpp" }
   --    defines {"NEA=16", "MAX_DAYS=100", "NA=9"}
   --    flags {"ExtraWarnings"}

   --    configuration "Release"
   --       flags { "OptimizeSize", "OptimizeSpeed" }

   --    configuration "Debug"
   --       defines { "DEBUG" }
   --       flags { "Symbols" }

   --    configuration "Optimized"
   --       flags { "OptimizeSize", "OptimizeSpeed", "EnableSSE", "EnableSSE2" }
   ----------------------------------------------------------------------------
   -- project "spatial_covid"
   --    kind "ConsoleApp"
   --    language "C++"
   --    files { "../src/spatial_covid.cpp" }
   --    defines {"NEA=16", "MAX_DAYS=100", "NA=9"}
   --    flags {"ExtraWarnings"}

   --    configuration "Release"
   --       flags { "OptimizeSize", "OptimizeSpeed" }

   --    configuration "Debug"
   --       defines { "DEBUG" }
   --       flags { "Symbols" }

   --    configuration "Optimized"
   --       flags { "OptimizeSize", "OptimizeSpeed", "EnableSSE", "EnableSSE2" }
   -- ----------------------------------------------------------------------------
   -- project "spatial_covid0d"
   --    kind "ConsoleApp"
   --    language "C++"
   --    files { "../src/spatial_covid0d.cpp" }
   --    defines {"NEA=16", "MAX_DAYS=100", "NA=9"}
   --    flags {"ExtraWarnings"}

   --    configuration "Release"
   --       flags { "OptimizeSize", "OptimizeSpeed" }

   --    configuration "Debug"
   --       defines { "DEBUG" }
   --       flags { "Symbols" }

   --    configuration "Optimized"
   --       flags { "OptimizeSize", "OptimizeSpeed", "EnableSSE", "EnableSSE2" }
   -- ----------------------------------------------------------------------------
   -- project "spatial_covid2d"
   --    kind "ConsoleApp"
   --    language "C++"
   --    files { "../src/spatial_covid2d.cpp" }
   --    defines {"NEA=16", "MAX_DAYS=100", "NA=9"}
   --    flags {"ExtraWarnings"}

   --    configuration "Release"
   --       flags { "OptimizeSize", "OptimizeSpeed" }

   --    configuration "Debug"
   --       defines { "DEBUG" }
   --       flags { "Symbols" }

   --    configuration "Optimized"
   --       flags { "OptimizeSize", "OptimizeSpeed", "EnableSSE", "EnableSSE2" }
   -- ----------------------------------------------------------------------------
  
