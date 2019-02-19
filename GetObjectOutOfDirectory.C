TObject* GetObjectOutOfCollection(TCollection* coll, TString objectName, TObject* dirofSearchedFile = 0x0);

TObject* GetObjectOutOfDirectory(TString fileName, TString objectName, TObject* dirofSearchedFile = 0x0) {
  TFile* f = new TFile(fileName.Data());
  return GetObjectOutOfDirectory(f,objectName);
}

TObject* GetObjectOutOfDirectory(TDirectoryFile* d, TString objectName, TObject* dirofSearchedFile = 0x0) {
  if (!d || d->IsZombie())
    return 0x0;
  
  TObject* obj = 0x0;
  d->GetObject(objectName.Data(),obj);
  if (!obj) {
    TList* list = d->GetListOfKeys();
    for (Int_t n=0;n<list->GetEntries();n++) {
      TObject* dirobj = 0x0;
      d->GetObject(list->At(n)->GetName(),dirobj);
      if (!dirobj)
        continue;
      
      if (dirobj->InheritsFrom("TDirectory")) 
        obj = GetObjectOutOfDirectory((TDirectoryFile*)dirobj, objectName, dirofSearchedFile);
      else if (dirobj->InheritsFrom("TCollection"))
        obj = GetObjectOutOfCollection((TCollection*)dirobj, objectName, dirofSearchedFile);
      
      if (obj)
        break;
    }
  }
  else {
    dirofSearchedFile = d;
  }  
  
  return obj;
}

TObject *GetObjectOutOfCollection(TCollection* coll, TString objectName, TObject* dirofSearchedFile) {
  TObject* obj = 0x0;
  obj = coll->FindObject(objectName.Data());
  
  if (!obj) {
    TIter it = TIter(coll);
    TObject* collobj = it();
    while (collobj) {
      if (collobj->InheritsFrom("TDirectory")) 
        obj = GetObjectOutOfDirectory((TDirectoryFile*)collobj, objectName, dirofSearchedFile);
      else if (collobj->InheritsFrom("TCollection"))
        obj = GetObjectOutOfCollection((TCollection*)collobj, objectName, dirofSearchedFile);
      
      if (obj)
        break;
      
      collobj = it();
    }
  }
  else {
    dirofSearchedFile = coll;
  }
  return obj;
}