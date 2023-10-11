import sys, os, tarfile, itertools, shutil
try:
  import yaml
except ImportError:
  print ("Trying to Install required module: requests\n")
  os.system('pip3 install PyYAML')
import yaml

# Cartesian product of lists
def cartesian_prod(list1 = [], bool_var1 = False, list2 = [], bool_var2 = False, list3 = [], bool_var3 = False):
  full_list = []
  if bool_var1 and bool_var2 and bool_var3:
    full_list = list(itertools.product(list1,list2, list3))
  elif bool_var1 and bool_var2:
    full_list = list(itertools.product(list1,list2))
  elif bool_var1 and bool_var3:
    full_list = list(itertools.product(list1,list3))
  elif bool_var2 and bool_var3:
    full_list = list(itertools.product(list2,list3))

  return full_list

def moveAllFilesinDir(srcDir, dstDir):
  file_names = os.listdir(srcDir)

  # Check if both the are directories
  if os.path.isdir(srcDir) and os.path.isdir(dstDir):
    for file_name in file_names:
      shutil.move(os.path.join(srcDir, file_name), dstDir)
  else:
      print("srcDir & dstDir should be Directories")
  
  if os.path.isdir(srcDir):
     os.rmdir(srcDir) 


if __name__ == '__main__':
############## inputs #####################
  with open('rendered_wano.yml') as file:
          wano_file = yaml.full_load(file)
###########################################

  outdict = {}
  outdict["iter"] = []
  var_float = []
  var_int = []
  var_str = []
  keys = []

  if wano_file["Float"]:
      var_begin = wano_file["VarF-begin"]
      var_end = wano_file["VarF-end"]
      var_npoints = wano_file["N-points"]
      keys.append("Float")
      step = (var_end-var_begin)/var_npoints
      for i in range(0,var_npoints):
          z_0 = var_begin + step*i
          z_0 = round(z_0,3)
          var_float.append(float(z_0))    
      
  if wano_file["Int"]:
      var_begin = wano_file["VarI-begin"]
      var_end = wano_file["VarI-end"]
      step = wano_file["Step"]
      if step == 0:
          print( "N_points must fit in the interval" )
          sys.exit()
          
      for z_0 in range(var_begin,var_end, step):
          #z_0 = var_begin + int(step*i)
          print(z_0)
          var_int.append(z_0)
      keys.append("Int")

  if wano_file["Structures"]:
    if os.path.isfile("Structures.tar"):
      print ("File exist")
      file = tarfile.open("Structures.tar",'r')
      file.extractall('Structures')
      #[print(os.path.isdir('Structures/'+ jj)) for jj in var_str] 


      srcDir = 'Structures/Molecules'
      dstDir = 'Structures'
      moveAllFilesinDir(srcDir, dstDir)
      var_str = os.listdir(dstDir)

      keys.append("Structures")  
      #outdict["iter"] = file.getnames()
      #print(outdict)
  else:
    print ("File not exist")
    tar = tarfile.open("Structures.tar", "w")
    tar.close()
    var_str = []   

  # Catersian product
  true_list = [wano_file["Float"], wano_file["Int"], wano_file["Structures"]]
  ii = 0
  for jj in true_list:
    if jj:
      ii += 1
  
  if ii > 1:
    input_list = cartesian_prod(var_float,true_list[0], var_int, true_list[1], var_str, true_list[2])
  elif true_list[0]:
    input_list = var_float
  elif true_list[1]:
    input_list = keys
  elif true_list[2]:
    input_list = var_str  

  
  #input_list = list(itertools.product(var_float,var_str))
  #print(input_list)
  
  # input_list = [list(x) for x in input_list]

  input_dict = {}  
  temp_dict = {}
  keys1 = list(range(len(input_list)))
  outdict["iter"] = keys1
  
  for i in keys1:
      temp_dict = {i: dict(zip(keys,input_list[i]))}
      input_dict.update(temp_dict)
  

  with open("input_dict.yml",'w') as out:
    yaml.dump(input_dict, out, default_flow_style = False)

  with open("output_dict.yml",'w') as out:
    yaml.dump(outdict, out,default_flow_style = False)