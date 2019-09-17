class TF1;

class PieceWisePoly {
public:
  PieceWisePoly(Int_t parts, Double_t* cutxvalues, Int_t* polys, Double_t xmin = 0, Double_t xmax = 1, Double_t* params = 0x0, Int_t smooth = 2);
  ~PieceWisePoly();
  double operator() (double* x, double* p = 0x0);
  double Evaluate(double* x, double* p = 0x0);
  double Eval (double x, double* p = 0x0);
  void SetParam(Double_t* params);
  Double_t SumUp(Double_t constant, Int_t start, Int_t end, Int_t derivative, Int_t startn);
  Int_t GetNOfParam() const {return nOfFreeParams;};
  TF1* GetFunction() {return piecewisepolynom;};
  TF1* GetPartFunction(Int_t i);
  Int_t GetNParts() {return nParts;};
  Double_t* GetCuts() {return cuts;};
  Int_t* GetNParameters() {return polyParameters;};
  
private:
  Int_t doSmoothing;
  Int_t nParts;
  Double_t* cuts;
  Int_t nOfFreeParams;
  Int_t* polyParameters;
  TF1* piecewisepolynom;
};


