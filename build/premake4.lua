
function queryTerminal(command)
   local success, handle = pcall(io.popen, command)
   if not success then 
       return ""
   end

   result = handle:read("*a")
   handle:close()
   result = string.gsub(result, "\n$", "") -- remove trailing whitespace
   return result
end

function getPythonPath()
   return queryTerminal("python3 -m pybind11 --includes")
end

function getPythonExtention()   
   return queryTerminal("python3-config --extension-suffix")
end

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
      buildoptions {  "-fPIC" }

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
      
   project "python_model_binding"
      pythonIncludes = getPythonPath()
      buildoptions {  pythonIncludes, "-fPIC" }
      kind "SharedLib"
      language "C++"
      targetname("cmodels") -- this name must match the module name in the macro PYBIND11_MODULE(group6_pybind_test, m)
      targetextension( getPythonExtention() )
      files ({
            "../src/python_interface.cpp"
         }
      )
      targetprefix("")
   

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
  
