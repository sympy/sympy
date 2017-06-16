
Map2[func_,lst1_,lst2_] :=
  Module[{i},
    ReapList[Do[Sow[func[lst1[[i]],lst2[[i]]]],{i,Length[lst1]}]]]
