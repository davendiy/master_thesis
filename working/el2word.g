# ElementAsWordGeneratorsPointGroup := function(S,g)
#   local gen,F,hom;
#   gen := GeneratorsOfGroup(S);
#   if Size(gen) = 1 then F := FreeGroup("a"); fi;
#   if Size(gen) = 2 then F := FreeGroup("a","b"); fi;
#   if Size(gen) = 3 then F := FreeGroup("a","b","c"); fi;
#   if Size(gen) = 4 then F := FreeGroup("a","b","c","d"); fi;
#
#   hom := GroupHomomorphismByImages(F,S,GeneratorsOfGroup(F),GeneratorsOfGroup(S));
#   if g in S then return PreImagesRepresentative(hom,g);
#   else Print("Element not in the group"); return false; fi;
# end;;





ElementAsWordGeneratorsPointGroup := function(S, g, F)
  local gen, hom;
  gen := MinimalGeneratingSet(S);
  hom := GroupHomomorphismByImages(F,S, MinimalGeneratingSet(F), MinimalGeneratingSet(S));
  if g in S then return PreImagesRepresentative(hom, g);
  else Print("Element not in the group"); return false; fi;
end;;


ElementAsWordGenerators:=function(G,g)
  local gen,S,genS,F_S,F_G,F_T,homS,
        homG,homSG,homT,homTG,s,word,wordS,
        wordT,n1,n2,g_linear,g_transl,e,f;

  if not (g in G) then
    Print("Element not in the group: ");
    return false;
  fi;

  gen := MinimalGeneratingSet(G);
  S := PointGroup(G); genS := MinimalGeneratingSet(S);
  if Size(genS) = 0 then
    n1 := g[3][1];
    n2 := g[3][2];
    F_T := FreeGroup("e","f"); e := F_T.1; f := F_T.2;
    word:=e^(n1)*f^(n2);
    return word;
  fi;

  if Size(genS)=1 then
    F_S:=FreeGroup("a");
    F_G:=FreeGroup("a","e","f");
  fi;

  if Size(genS)=2 then
    F_S:=FreeGroup("a","b");
    F_G:=FreeGroup("a","b","e","f");
  fi;

  if Size(genS)=3 then
    F_S:=FreeGroup("a","b","c");
    F_G:=FreeGroup("a","b","c","e","f");
  fi;

  if Size(genS)=4 then
    F_S:=FreeGroup("a","b","c","d");
    F_G:=FreeGroup("a","b","c","d","e","f");
  fi;

  homS:=GroupHomomorphismByImages(F_S,S,GeneratorsOfGroup(F_S),GeneratorsOfGroup(S));
  homG:=GroupHomomorphismByImages(F_G,G,GeneratorsOfGroup(F_G),GeneratorsOfGroup(G));
  homSG:=GroupHomomorphismByImages(F_S,F_G,GeneratorsOfGroup(F_S),GeneratorsOfGroup(F_G){[1..Size(GeneratorsOfGroup(F_S))]});

  F_T:=FreeGroup("e","f"); e:=F_T.1; f:=F_T.2;
  homT:=GroupHomomorphismByImages(F_T,G,GeneratorsOfGroup(F_T),gen{[Size(GeneratorsOfGroup(S))+1..Size(GeneratorsOfGroup(S))+2]});
  homTG:=GroupHomomorphismByImages(F_T,F_G,GeneratorsOfGroup(F_T),GeneratorsOfGroup(F_G){[Size(GeneratorsOfGroup(F_S))+1..Size(GeneratorsOfGroup(F_S))+2]});

  s:=g{[1,2]}{[1,2]};
  wordS:=PreImagesRepresentative(homS,s);
  g_linear:=(wordS^homSG)^homG;
  n1:=(g*g_linear^-1)[3][1];
  n2:=(g*g_linear^-1)[3][2];
  g_transl:=gen[Size(genS)+1]^(n1)*gen[Size(genS)+2]^(n2);
    wordT:=e^(n1)*f^(n2);
  word:=wordS^homSG*wordT^homTG;

  if (g_linear * g_transl <> g) then Print("Error"); fi;
  return word;
end;;
