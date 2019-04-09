#include "FJ_includes.h"
#include "AliFJWrapper.h"

void Test() {
//   AliFJWrapper* wrap = new AliFJWrapper();
  TString name = "wrap1";
  TString title = "wrap12";
  AliFJWrapper* wrap = new AliFJWrapper(name.Data(),title.Data());
  printf("%s",wrap->GetName());
}