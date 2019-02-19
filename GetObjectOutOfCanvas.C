#include "TObject.h"
#include "TCanvas.h"
#include "TList.h"
#include <iostream>
using namespace std;

TObject *GetObjectOutOfListIterative(TList *list, const char *name);

TObject* GetObjectOutOfCanvas(TPad *p, const char *name) { 
  return GetObjectOutOfListIterative(p->GetListOfPrimitives(), name);
}

TObject *GetObjectOutOfListIterative(TList *list, const char *name) {
  TObject *obj = 0x0;
  obj = list->FindObject(name);
  if (obj) 
   return obj;
  
  Int_t ii = 0;
  TObject *iterobj = list->At(ii);
  while(iterobj) {
    TString cname(iterobj->ClassName());
    if (!cname.CompareTo("TPad") || !cname.CompareTo("TCanvas")) {
      obj = GetObjectOutOfCanvas((TPad*)iterobj, name);
      if (obj) 
        return obj;
    }
    if (!cname.CompareTo("TList")) {
      obj = GetObjectOutOfListIterative((TList*)iterobj, name);
      if (obj)
        return obj;
    }
    ii++;
    iterobj = list->At(ii);
  }
  return 0x0;
}