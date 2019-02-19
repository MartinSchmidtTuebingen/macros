void Test() {
  gSystem->AddIncludePath("-$FASTJET_ROOT/include");
  gSystem->Load("libfastjet");
  AliFJWrapper* wrap = new AliFJWrapper();
}