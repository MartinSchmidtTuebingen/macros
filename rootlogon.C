using namespace std;
#include "TStyle.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TObject.h"
#include "TString.h"
#include "TLine.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TGraph.h"
#include "iostream"
#include "TF1.h"
#include "TAxis.h"
#include "TMultiGraph.h"

Int_t CheckLoadLibrary(const char* library);

void rootlogon()
{
  TString aliceRoot = gSystem->Getenv("ALICE_ROOT");
  Bool_t alirootFound = kTRUE;
  if (aliceRoot == "") {
    alirootFound = kFALSE;
    printf("Aliroot not loaded. Only loading root libraries...\n");
  }
  
  //Old AliRoot gSystem->AddDynamicPath("/home/administrator/Alice/aliroot/master/build/lib/tgt_linuxx8664gcc");
  if (alirootFound)
    gSystem->AddDynamicPath("$ALICE_ROOT/lib");
  
  //set inlcude paths
  if (alirootFound) {
    gSystem->AddIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/ -I$ALICE_ROOT/include -I$ALICE_ROOT/STEER  -I$ALICE_ROOT/ANALYSIS -I$ALICE_ROOT/PWGCF -I$ALICE_ROOT/PWGDQ -I$ALICE_ROOT/PWGGA -I$ALICE_ROOT/PWGHF -I$ALICE_ROOT/PWGJE -I$ALICE_ROOT/PWGLF -I$ALICE_ROOT/PWGPP -I$ALICE_ROOT/PWGUD -I$ALICE_PHYSICS/CORRFW -I$ALICE_ROOT/TPC -I$ALICE_ROOT/TPC/Base -I$ALICE_ROOT/TPC/Sim -I$ALICE_ROOT/TPC/Rec -I$ALICE_ROOT/TPC/Calib -I$ALICE_ROOT/TPC/Util -I$ALICE_ROOT/TRD -I$ALICE_ROOT/PWG/muon -I$ALICE_ROOT/JETAN -I$ALICE_ROOT/ANALYSIS/Tender -I$ALICE_ROOT/EMCAL");
  }
  else {
    gSystem->AddIncludePath("-I$ROOTSYS/include");
  }
  
  CheckLoadLibrary("libCore");
  CheckLoadLibrary("libPhysics");
  CheckLoadLibrary("libMinuit");
  CheckLoadLibrary("libGui");
  CheckLoadLibrary("libXMLParser");
  
  CheckLoadLibrary("libGeom");
  CheckLoadLibrary("libVMC");

  CheckLoadLibrary("libNet");
  CheckLoadLibrary("libTree");
  CheckLoadLibrary("libProof");

  if (alirootFound) {
    CheckLoadLibrary("libSTEERBase");
    CheckLoadLibrary("libESD");
    CheckLoadLibrary("libCDB");
    CheckLoadLibrary("libRAWDatabase");
    CheckLoadLibrary("libRAWDatarec");
    CheckLoadLibrary("libANALYSIS");
    CheckLoadLibrary("libSTEER");

    CheckLoadLibrary("libSTAT");

    CheckLoadLibrary("libAOD");
    CheckLoadLibrary("libOADB");
    CheckLoadLibrary("libANALYSISalice");
    CheckLoadLibrary("libCORRFW");
    CheckLoadLibrary("libTPCbase");
  }
  
  CheckLoadLibrary("libCint");//To fix bug in new root version

  
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  mystyle();

  if (alirootFound) {
    gROOT->ProcessLine(".include $ALICE_ROOT/include");
    gROOT->ProcessLine(".include $ALICE_PHYSICS/include");
    gROOT->ProcessLine(".include $ALICE_ROOT/TRD");
  }
  printf("rootlogon.C executed\n");
}

Int_t CheckLoadLibrary(const char* library)
{
  // checks if a library is already loaded, if not loads the library

  if (strlen(gSystem->GetLibraries(Form("%s", library), "", kFALSE)) > 0)
    return 1;

  return gSystem->Load(library);
}

void  mystyle()
{
    // to compile the lib do once
    // .L mystyle.C++
    // if lib has been compiled load the lib
    // add this lines to your macro
    // gSystem->Load("/misc/pachmay/bin/mystyle_C.so");
    // mystyle();


    cout<<"LOADING PERSONAL STYLE SETTING!"<<endl;

    Float_t titlesize=0.05;
//    Float_t titleoffset=1.25;
//    Float_t titleoffsety=1.5;
    Float_t titleoffset=0.9;
    Float_t titleoffsety=0.9;
    Float_t labelsize=0.05;


    // defining canvas options ------------------------
    gStyle->SetCanvasColor(10);
    gStyle->SetPadColor(10);
    
    // Set Ticks on all sides (alternatively: TCanvas->SetTickx() / ->SetTicky())
    // 2 instead of 1 means that there are also labels on the other side (e.g. gStyle->SetPadTickY(2))
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    
    // Set grid in x and y
    gStyle->SetPadGridX(1);
    gStyle->SetPadGridY(1);
/*
    // defining size of pad ---------------------------
    gStyle->SetPadTopMargin(0.15);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetPadLeftMargin(0.15);
*/
    gStyle->SetPadBottomMargin(0.14);
    gStyle->SetPadBorderMode(0);

    // defining outfit of axis ------------------------
    // offset of titles
    gStyle->SetTitleXOffset(titleoffset);
    gStyle->SetTitleYOffset(titleoffsety);
    // size of titles
    gStyle->SetTitleSize(titlesize,"X");
    gStyle->SetTitleSize(titlesize,"Y");
    gStyle->SetTitleSize(titlesize,"Z");
    // size of axis labels
    gStyle->SetLabelSize(labelsize,"X");
    gStyle->SetLabelSize(labelsize,"Y");
    gStyle->SetLabelSize(labelsize,"Z");

    // switch options of the histogram ----------------
    //gStyle->SetOptStat(0);
    gStyle->SetOptStat(111111);   // switch on(1)/off(0) stat box; 111111: Also print under- and overflow
    //gStyle->SetOptStat("emr"); //Only box with entries, mean and RMS
    gStyle->SetOptFit(1);         // show fit results
    gStyle->SetErrorX(0.5);          // length of error bar in x (take 0.5 for error equal bin size!)
    gStyle->SetTitleH(0.08);       // make title larger
    //gStyle->SetTitleW(0.8);       // make title larger
    gStyle->SetTitleStyle(0);      // make title box transparent
    gStyle->SetTitleBorderSize(0); // make border of box disapear
    
//     gROOT->LoadMacro("/home/administrator/aliceworkinprogress.C");
}

void setOPT_canvas(TCanvas *canvas)
{

    gStyle->SetLineStyleString(22,"80 18 12 18 12 12");
    gStyle->SetPadColor(10);

    gStyle->SetPalette(1);
    gStyle->SetOptDate(0);

    canvas->SetLeftMargin(0.2);
    canvas->SetRightMargin(0.15);
    canvas->SetTopMargin(0.15);
    canvas->SetBottomMargin(0.17);



    canvas->SetFillColor(10);
    canvas->SetBorderMode(0);
    canvas->SetBorderSize(0);
    canvas->SetTickx();
    canvas->SetTicky();
    canvas->SetFrameFillColor(0);
    canvas->SetFrameLineWidth(2);
    canvas->SetFrameBorderMode(0);
    canvas->SetFrameBorderSize(0);


}
void setDefaultOptions_hist(TH1F *hist, Bool_t clearTitle = kTRUE)
{
  hist->GetXaxis()->SetTitleOffset(0.9);
  hist->SetLineWidth(1);
  if (clearTitle) hist->SetTitle("");
}

void setOPT_hists(TH1F *hist, TString xAxisTitle=" ", TString yAxisTitle=" ", Int_t Ndevision=510, Int_t marker_style =20, Float_t marker_size =1.3, Int_t color=1)
{
    gStyle->SetPadRightMargin(0.05);
    hist->SetTitle("");
    hist->GetXaxis()->SetTitle(xAxisTitle);
    hist->GetYaxis()->SetTitle(yAxisTitle);
    hist->GetXaxis()->SetTitleSize(0.06);
    hist->GetYaxis()->SetTitleSize(0.06);
    hist->GetXaxis()->SetTitleFont(42);
    //hist->GetXaxis()->SetDecimals(kTRUE);
    //hist->GetYaxis()->SetDecimals(kTRUE);
    hist->GetYaxis()->SetTitleFont(42);
    hist->GetXaxis()->SetNdivisions(Ndevision);
    hist->GetYaxis()->SetTitleOffset(1.2);
    hist->GetXaxis()->SetTitleOffset(1.2);
//    hist->GetXaxis()->CenterTitle();
//    hist->GetYaxis()->CenterTitle();
    hist->GetXaxis()->SetLabelFont(42);
    hist->GetYaxis()->SetLabelFont(42);
    hist->GetXaxis()->SetLabelSize(0.05);
    hist->GetYaxis()->SetLabelSize(0.05);
    hist->SetMarkerStyle(marker_style);
    hist->SetMarkerSize(marker_size);
    hist->SetMarkerColor(color);
    hist->SetLineColor(color);
    hist->SetLineWidth(2);
}

void setOPT_hists2D(TH2F *hist2, TString xAxisTitle=" ", TString yAxisTitle=" ", Int_t Ndevision=510)
{
    gStyle->SetPadRightMargin(0.15);
    hist2->SetTitle("");
    hist2->GetXaxis()->SetTitle(xAxisTitle);
    hist2->GetYaxis()->SetTitle(yAxisTitle);
    hist2->GetXaxis()->SetTitleSize(0.06);
    hist2->GetYaxis()->SetTitleSize(0.06);
    hist2->GetZaxis()->SetTitleSize(0.06);
    hist2->GetXaxis()->SetTitleFont(42);
    hist2->GetYaxis()->SetTitleFont(42);
    hist2->GetZaxis()->SetTitleFont(42);
    hist2->GetXaxis()->SetNdivisions(Ndevision);
    hist2->GetYaxis()->SetTitleOffset(1.2);
    hist2->GetXaxis()->SetTitleOffset(1.2);
//    hist2->GetXaxis()->CenterTitle();
//    hist2->GetYaxis()->CenterTitle();
    hist2->GetXaxis()->SetLabelFont(42);
    hist2->GetYaxis()->SetLabelFont(42);
    hist2->GetZaxis()->SetLabelFont(42);
    hist2->GetXaxis()->SetLabelSize(0.05);
    hist2->GetYaxis()->SetLabelSize(0.05);
    hist2->GetZaxis()->SetLabelSize(0.05);
    hist2->GetXaxis()->SetDecimals(kTRUE);
    hist2->GetYaxis()->SetDecimals(kTRUE);
    hist2->GetZaxis()->SetDecimals(kTRUE);
}


void setOPT_histsF(TF1 *hist, TString xAxisTitle=" ", TString yAxisTitle=" ", Int_t Ndevision=510, Int_t marker_style =20, Float_t marker_size =1.3, Int_t color=1)
{
    gStyle->SetPadRightMargin(0.05);
    hist->GetXaxis()->SetDecimals(kTRUE);
    hist->GetXaxis()->SetTitle(xAxisTitle);
    hist->GetYaxis()->SetTitle(yAxisTitle);
    hist->GetXaxis()->SetTitleSize(0.06);
    hist->GetYaxis()->SetTitleSize(0.06);
    hist->GetXaxis()->SetTitleFont(42);
    hist->GetYaxis()->SetTitleFont(42);
    hist->GetXaxis()->SetNdivisions(Ndevision);
    hist->GetYaxis()->SetTitleOffset(1.2);
    hist->GetXaxis()->SetTitleOffset(1.2);
//    hist->GetXaxis()->CenterTitle();
//    hist->GetYaxis()->CenterTitle();
    hist->GetXaxis()->SetLabelFont(42);
    hist->GetYaxis()->SetLabelFont(42);
    hist->GetXaxis()->SetLabelSize(0.05);
    hist->GetYaxis()->SetLabelSize(0.05);
    hist->SetMarkerStyle(marker_style);
    hist->SetMarkerSize(marker_size);
    hist->SetMarkerColor(color);
    hist->SetLineColor(color);
    hist->SetLineWidth(2);
}

void setOPT_graph(TGraph *hist, TString xAxisTitle=" ", TString yAxisTitle=" ", Int_t Ndevision=510, Int_t marker_style =21, Float_t marker_size =1.3, Int_t color=1)
{
    gStyle->SetPadRightMargin(0.05);
    hist->GetYaxis()->SetDecimals(kTRUE);
    hist->GetXaxis()->SetDecimals(kTRUE);
    hist->GetXaxis()->SetTitle(xAxisTitle);
    hist->GetYaxis()->SetTitle(yAxisTitle);
    hist->GetXaxis()->SetTitleSize(0.06);
    hist->GetYaxis()->SetTitleSize(0.06);
    hist->GetYaxis()->SetTitleOffset(1.2);
    hist->GetXaxis()->SetTitleOffset(1.2);
    hist->GetXaxis()->SetTitleFont(42);
    hist->GetYaxis()->SetTitleFont(42);
    hist->GetXaxis()->SetNdivisions(Ndevision);
 
    
    hist->GetXaxis()->SetLabelFont(42);
    hist->GetYaxis()->SetLabelFont(42);
    hist->GetXaxis()->SetLabelSize(0.05);
    hist->GetYaxis()->SetLabelSize(0.05);
    hist->SetMarkerStyle(marker_style);
    hist->SetMarkerSize(marker_size);
    hist->SetMarkerColor(color);
    hist->SetLineColor(color);
}

void setOPT_multigraph(TMultiGraph *hist, TString xAxisTitle=" ", TString yAxisTitle=" ", TString title=" ")
{
    gStyle->SetPadRightMargin(0.05);
    hist->SetTitle(title);
    hist->GetYaxis()->SetDecimals(kTRUE);
    hist->GetXaxis()->SetDecimals(kTRUE);
    hist->GetXaxis()->SetTitle(xAxisTitle);
    hist->GetYaxis()->SetTitle(yAxisTitle);
    hist->GetXaxis()->SetTitleSize(0.06);
    hist->GetYaxis()->SetTitleSize(0.06);
    hist->GetYaxis()->SetTitleOffset(1.2);
    hist->GetXaxis()->SetTitleOffset(1.2);
    hist->GetXaxis()->SetTitleFont(42);
    hist->GetYaxis()->SetTitleFont(42);
    hist->GetXaxis()->SetLabelFont(42);
    hist->GetYaxis()->SetLabelFont(42);
    hist->GetXaxis()->SetLabelSize(0.05);
    hist->GetYaxis()->SetLabelSize(0.05);
   
}

TLegend* plotLegend(TString pos="right_top",TString Title="No Title",Float_t scaleX=0.9,Float_t scaleY=0.9,Float_t offsetX=0.0,Float_t offsetY=0.0,TString Comment="",Int_t commentcolor=1)
{
    // pos places the Legend and can be
    // right_top, mid_top, left_top, right_bottom,mid_bottom,left_bottom
    // Title gives you the Title of the legend
    // scaleX,scaleY defines the size (=1 means 1/3 pad with,1/3 pad height (according to pos))
    // offsetX , offsetY (NDC coordinates of pad -> 0.1) shift by 10% of pad with) shifts the legend by the value in x and y
    // comment is optional text below title before legend entris will start
    // commentcolor defines the text color of the comment

    Float_t left  =gPad->GetLeftMargin()*1.15;
    Float_t right =1-gPad->GetRightMargin();
    Float_t top   =1-gPad->GetTopMargin();
    Float_t bottom=gPad->GetBottomMargin()*1.15;
    Float_t mid   =gPad->GetLeftMargin() + (1-(gPad->GetRightMargin()+gPad->GetLeftMargin()))/2;
    Float_t width =(right-left)/2;
    Float_t heith =(top-bottom)/2;
    TLegend* legend;
    TLine* dummy=new TLine();
    dummy->SetLineColor(10);

    if(pos.CompareTo("left_top")==0)    legend=new TLegend(left+offsetX,top+offsetY-(scaleY*heith),left+offsetX+(scaleX*width),top+offsetY,Title.Data());                           // left top
    if(pos.CompareTo("mid_top")==0)     legend=new TLegend(mid+offsetX-((scaleX*width)/2),top+offsetY-(scaleY*heith),mid+offsetX+((scaleX*width)/2),top+offsetY,Title.Data());      // mid up
    if(pos.CompareTo("right_top")==0)   legend=new TLegend(right+offsetX-(scaleX*width),top+offsetY-(scaleY*heith),right+offsetX,top+offsetY,Title.Data());                         // right top
    if(pos.CompareTo("left_bottom")==0) legend=new TLegend(left+offsetX,bottom+offsetY,left+offsetX+(scaleX*width),bottom+offsetY+(scaleY*heith),Title.Data());                     // left bottom
    if(pos.CompareTo("mid_bottom")==0)  legend=new TLegend(mid+offsetX-((scaleX*width)/2),bottom+offsetY,mid+offsetX+((scaleX*width)/2),bottom+offsetY+(scaleY*heith),Title.Data());// mid down
    if(pos.CompareTo("right_bottom")==0)legend=new TLegend(right+offsetX-(scaleX*width),bottom+offsetY,right+offsetX,bottom+offsetY+(scaleY*heith),Title.Data());                   // right bottom

    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetMargin(0.15);



    if(Comment.CompareTo("")!=0)
    {   // comment below header
        TLegendEntry* entry=legend->AddEntry(dummy,Comment.Data());
        entry->SetTextColor(commentcolor);
        entry->SetTextFont(62);
        entry->SetOption("l");
    }
    return legend;
}
void setLegendEntry(TLegend* legend,TObject* object,char* label,Int_t col,TString opt,Float_t size=0.05)
{
    // add entry for object with label and color and option to legend
    // legend is thge pointer to the legend where you want to add the entry
    // object ist the histogram or graph which will be added
    // label is the text which will appear in the entry
    // col is th color of the text
    // opt is the option ("pL"(marker and line), "p"(marker),"l"(line))
    TLegendEntry* entry=legend->AddEntry((TObject*)object,label);
    entry->SetTextColor(col);
    entry->SetTextSize(size);
    entry->SetOption(opt.Data());
}
TLatex* plotTopLegend(char* label,Float_t x=-1,Float_t y=-1,Float_t size=0.05,Int_t color=1)
{
    // coordinates in NDC!
    // plots the string label in position x and y in NDC coordinates
    // size is the text size
    // color is the text color

    if(x<0||y<0)
    {   // defaults
        x=gPad->GetLeftMargin()*1.15;
        y=(1-gPad->GetTopMargin())*1.02;
    }
    TLatex* text=new TLatex(x,y,label);
    text->SetTextSize(size);
    text->SetNDC();
    text->SetTextColor(color);
    text->Draw();
    return text;
}
// void saveToPdf(TCanvas* canvas,Char_t* path,Char_t* filename,Bool_t separate=kFALSE)
// {
//     // saves a canvas to pdf file by storing it as ps and calling
//     // ps2pdf to convert it to pdf. The ps file will be removed.
//     // canvas is the pointer of the canvas which you want to save
//     // path is the path to the directory where you want to store the pdf
//     // filename is the name of the pdf which wou create
//     // separate=kTRUE will loop over all pads and save the pictures single (path/filename_i)
//     // separate=kFALSE as the canvas in one pdf.
//     TString File    =filename;
//     TString Path    =path;
//     Path            =Path + "/";
//     TString command ="";
//     TString totalin ="";
//     TString totalout="";

//     File.ReplaceAll(".ps","");
//     File.ReplaceAll(".pdf","");

//     if(!separate)
//     {
//         totalin = Path + File + ".ps";
//         totalout= Path + File + ".pdf";
//         canvas->SaveAs(totalin.Data());
//         command= "ps2pdf " + totalin + " " + totalout;
//         gSystem->Exec(command.Data());
//         command= "rm " + totalin;
//         gSystem->Exec(command.Data());
//     }else{
//         TString subpad="";
//         Int_t npads= canvas->GetListOfPrimitives()->GetSize();
//         for(Int_t i=1;i<=npads;i++)
//         {
//             subpad=i;
//             totalin = Path + File + "_" + subpad + ".ps";
//             totalout= Path + File + "_" + subpad + ".pdf";
//             canvas->cd(i);
//             gPad->SaveAs(totalin.Data());
//             command= "ps2pdf " + totalin + " " + totalout;
//             gSystem->Exec(command.Data());
//             command= "rm " + totalin;
//             gSystem->Exec(command.Data());
//         }
//     }
// }
// void saveToPng(TCanvas* canvas,Char_t* path,Char_t* filename,Bool_t separate=kFALSE)
// {
//     // saves a canvas to png file by storing it as gif and calling
//     // convert to convert it to png. The gif file will be removed.
//     // canvas is the pointer of the canvas which you want to save
//     // path is the path to the directory where you want to store the png
//     // filename is the name of the png which wou create
//     // separate=kTRUE will loop over all pads and save the pictures single (path/filename_i)
//     // separate=kFALSE as the canvas in one png.
//     TString File    =filename;
//     TString Path    =path;
//     Path            =Path + "/";
//     TString command ="";
//     TString totalin ="";
//     TString totalout="";


//     File.ReplaceAll(".ps","");
//     if(!separate)
//     {
//         totalin = Path + File + ".ps";
//         totalout= Path + File + ".png";
//         canvas->SaveAs(totalin.Data());
//         command= "convert " + totalin + " " + totalout;
//         gSystem->Exec(command.Data());
//         command= "rm " + totalin;
//         gSystem->Exec(command.Data());
//     }else{
//         TString subpad="";
//         Int_t npads= canvas->GetListOfPrimitives()->GetSize();
//         for(Int_t i=1;i<=npads;i++)
//         {
//             subpad=i;
//             totalin = Path + File + "_" + subpad + ".ps";
//             totalout= Path + File + "_" + subpad + ".png";
//             canvas->cd(i);
//             gPad->SaveAs(totalin.Data());
//             command= "convert " + totalin + " " + totalout;
//             gSystem->Exec(command.Data());
//             command= "rm " + totalin;
//             gSystem->Exec(command.Data());
//         }
//     }
// }

// void saveToC(TCanvas* canvas,Char_t* path,Char_t* filename,Bool_t separate=kFALSE)
// {
//     // saves a canvas to png file by storing it as gif and calling
//     // convert to convert it to png. The gif file will be removed.
//     // canvas is the pointer of the canvas which you want to save
//     // path is the path to the directory where you want to store the png
//     // filename is the name of the png which wou create
//     // separate=kTRUE will loop over all pads and save the pictures single (path/filename_i)
//     // separate=kFALSE as the canvas in one png.
//     TString File    =filename;
//     TString Path    =path;
//     Path            =Path + "/";
//     TString command ="";
//     TString totalin ="";
//     TString totalout="";

//     File.ReplaceAll(".C","");
//     if(!separate)
//     {
//         totalin = Path + File + ".C";
//         canvas->SaveAs(totalin.Data());
//     }else{
//         TString subpad="";
//         Int_t npads= canvas->GetListOfPrimitives()->GetSize();
//         for(Int_t i=1;i<=npads;i++)
//         {
//             subpad=i;
//             totalin = Path + File + "_" + subpad + ".C";
//             canvas->cd(i);
//             gPad->SaveAs(totalin.Data());
//         }
//     }
// }

// void saveGifToPdf(TCanvas* canvas,Char_t* path,Char_t* filename,Bool_t separate=kFALSE)
// {
//     // saves a canvas to pdf file by storing it as gif and calling
//     // convert to convert it to pdf. The gif file will be removed.
//     // canvas is the pointer of the canvas which you want to save
//     // path is the path to the directory where you want to store the pdf
//     // filename is the name of the pdf which wou create
//     // separate=kTRUE will loop over all pads and save the pictures single (path/filename_i)
//     // separate=kFALSE as the canvas in one pdf.
//     TString File    =filename;
//     TString Path    =path;
//     Path            =Path + "/";
//     TString command ="";
//     TString totalin ="";
//     TString totalout="";

//     File.ReplaceAll(".pdf","");
//     if(!separate)
//     {
//         totalin = Path + File + ".gif";
//         totalout= Path + File + ".pdf";
//         canvas->SaveAs(totalin.Data());
//         command= "convert " + totalin + " " + totalout;
//         gSystem->Exec(command.Data());
//         command= "rm " + totalin;
//         gSystem->Exec(command.Data());
//     }else{
//         TString subpad="";
//         Int_t npads= canvas->GetListOfPrimitives()->GetSize();
//         for(Int_t i=1;i<=npads;i++)
//         {
//             subpad=i;
//             totalin = Path + File + "_" + subpad + ".gif";
//             totalout= Path + File + "_" + subpad + ".pdf";
//             canvas->cd(i);
//             gPad->SaveAs(totalin.Data());
//             command= "convert " + totalin + " " + totalout;
//             gSystem->Exec(command.Data());
//             command= "rm " + totalin;
//             gSystem->Exec(command.Data());
//         }
//     }
// }
Int_t getColor(Int_t i)
{
    // pre defined set of 10 colors
    // i is larger 9 it will start from the first color again
    Int_t colors[10]   ={2,4,6,8,38,46,14,1,30,43};
    return colors[i%10];
}
Int_t getMarker(Int_t i)
{
    // pre defined set of 10 marker styles
    // i is larger 9 it will start from the first marker style again
    Int_t style[10]   ={20,21,22,23,24,25,26,27,28,29};
    return style[i%10];
}
void setGraph(TGraph* graph,Int_t i,Int_t markerstyle=20,Int_t markercolor=1,Float_t markersize=1.,Int_t linecolor=1)
{
    // sets the graphics of graph (use full in a loop)
    // graph is the pointer of the graph
    //
    // i switches to modes:
    // (if i>0) i the number of the color of the marker and the markert style (getColor(),getMarker())
    // (if i<0) marker color and style will be  markercolor, markerstyle
    //
    // markersize, and linecolor defines the size of the marker and the line color
    // can be use in two ways:
    // values >=0 the style and color will be set via getColor(i),getMarker(i)
    // values < 0 the style and color will be set -markercolor -markerstyle

    if(i<0)
    {
        graph->SetMarkerColor(markercolor);
        graph->SetMarkerSize(markersize);
        graph->SetMarkerStyle(markerstyle);
        graph->SetLineColor(linecolor);
    }else{
        if     (markercolor>=0)graph->SetMarkerColor(getColor(i));
        else if(markercolor<0) graph->SetMarkerColor(getColor(-markercolor));
        graph->SetMarkerSize(markersize);
        if     (markerstyle>=0)graph->SetMarkerStyle(getMarker(i));
        else if(markerstyle<0) graph->SetMarkerStyle(getMarker(-markerstyle));
        graph->SetLineColor(linecolor);
    }
}
void setHist(TH1* hist,Int_t i,Int_t markerstyle=20,Int_t markercolor=1,Float_t markersize=1.,Int_t linecolor=1)
{
    // sets the graphics of histogram (use full in a loop)
    // hist is the pointer of the hiistogram
    //
    // i switches to modes:
    // (if i>0) i the number of the color of the marker and the markert style (getColor(),getMarker())
    // (if i<0) marker color and style will be  markercolor, markerstyle
    //
    // markersize, and linecolor defines the size of the marker and the line color
    // can be use in two ways:
    // values >=0 the style and color will be set via getColor(i),getMarker(i)
    // values < 0 the style and color will be set -markercolor -markerstyle
    if(i==-99)
    {
        hist->SetMarkerColor(markercolor);
        hist->SetMarkerSize(markersize);
        hist->SetMarkerStyle(markerstyle);
        hist->SetLineColor(linecolor);
    }else{
        if     (markercolor>=0)hist->SetMarkerColor(getColor(i));
        else if(markercolor<0) hist->SetMarkerColor(getMarker(-markercolor));
        hist->SetMarkerSize(markersize);
        if     (markerstyle>=0)hist->SetMarkerStyle(getMarker(i));
        else if(markerstyle<0) hist->SetMarkerStyle(getMarker(-markerstyle));
        hist->SetLineColor(linecolor);
      }
}

TH2D* histDiff(TH2D* h1, TH2D* h2, TString drawOption = "surf1")
{
  if (!h1 || !h2)
    return 0x0;
  
  TH2D* hDiff = (TH2D*)h1->Clone(Form("hDiff_%s_minus_%s", h1->GetName(), h2->GetName()));
  hDiff->Add(h2, -1);
  
  TCanvas* c = new TCanvas();
  c->SetTitle(Form("%s minus %s", h1->GetName(), h2->GetName()));
  c->cd();
  
  hDiff->DrawCopy(drawOption.Data());
  
  return hDiff;
}

TH2D* histRatio(TH2D* h1, TH2D* h2, TString drawOption = "surf1")
{
  if (!h1 || !h2)
    return 0x0;
  
  TH2D* hRatio = (TH2D*)h1->Clone(Form("hRatio_%s_divided_%s", h1->GetName(), h2->GetName()));
  hRatio->Divide(h2);
  
  TCanvas* c = new TCanvas();
  c->SetTitle(Form("%s divided by %s", h1->GetName(), h2->GetName()));
  c->cd();
  
  hRatio->DrawCopy(drawOption.Data());
  
  return hRatio;
}

TH2D* histRatioInterpolate(TH2D* h1, TH2D* h2, TString drawOption = "surf1")
{
  if (!h1 || !h2)
    return 0x0;
  
  TH2D* hRatio = (TH2D*)h1->Clone(Form("hRatio_%s_divided_%s", h1->GetName(), h2->GetName()));
  
  for (Int_t binX = 1; binX <= hRatio->GetNbinsX(); binX++) {
    for (Int_t binY = 1; binY <= hRatio->GetNbinsY(); binY++) {
      Double_t binCenterX = h1->GetXaxis()->GetBinCenter(binX);
      Double_t binCenterY = h1->GetYaxis()->GetBinCenter(binY);
      
      hRatio->SetBinContent(binX, binY, h1->Interpolate(binCenterX, binCenterY) / h2->Interpolate(binCenterX, binCenterY));
    }
  }
  
  TCanvas* c = new TCanvas();
  c->SetTitle(Form("%s divided by %s", h1->GetName(), h2->GetName()));
  c->cd();
  
  hRatio->DrawCopy(drawOption.Data());
  
  return hRatio;
}


void compareFractionsUsingPiPlusMu(Bool_t plotIdentifiedSpectra)
{
  TCanvas* cFractions = new TCanvas("cFractions", "Particle fractions", 100, 10, 1200, 800);
  cFractions->SetGridx(1);
  cFractions->SetGridy(1);
  TH1F* hFractionPionsPlusMuons = new TH1F(*hFractionPions);
  hFractionPionsPlusMuons->SetName("hFractionPionsPlusMuons");
  hFractionPionsPlusMuons->Add(hFractionMuons);
  
  hFractionPionsPlusMuons->GetYaxis()->SetRangeUser(0.0, 1.0);
  hFractionPionsPlusMuons->GetXaxis()->SetMoreLogLabels(kTRUE);
  hFractionPionsPlusMuons->GetXaxis()->SetNoExponent(kTRUE);
  hFractionPionsPlusMuons->Draw("e p");
  
  TH1F* hFractionPionsPlusMuonsMC = 0x0;
  if (plotIdentifiedSpectra) {
    hFractionPionsPlusMuonsMC = new TH1F(*hFractionPionsMC);
    hFractionPionsPlusMuonsMC->SetName("hFractionPionsPlusMuonsMC");
    hFractionPionsPlusMuonsMC->Add(hFractionMuonsMC);
    hFractionPionsPlusMuonsMC->Draw("e p same");
  }
  
  hFractionKaons->Draw("e p same");
  if (plotIdentifiedSpectra) {
    hFractionKaonsMC->Draw("e p same");
  }
  
  hFractionProtons->Draw("e p same");
  if (plotIdentifiedSpectra) {
    hFractionProtonsMC->Draw("e p same");
  }
  
  hFractionElectrons->Draw("e p same");
  if (plotIdentifiedSpectra) {
    hFractionElectronsMC->Draw("e p same");
  }
  
  TLegend* legend = new TLegend(0.622126, 0.605932, 0.862069, 0.855932);    
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  if (plotIdentifiedSpectra)
    legend->SetNColumns(2);
  if (plotIdentifiedSpectra)
    legend->AddEntry((TObject*)0x0, "Fit", "");
  if (plotIdentifiedSpectra)
    legend->AddEntry((TObject*)0x0, "MC", "");
  legend->AddEntry(hFractionPionsPlusMuons, "#pi + #mu", "p");
  if (plotIdentifiedSpectra)
    legend->AddEntry(hFractionPionsPlusMuonsMC, "#pi + #mu", "p");
  legend->AddEntry(hFractionKaons, "K", "p");
  if (plotIdentifiedSpectra)
    legend->AddEntry(hFractionKaonsMC, "K", "p");
  legend->AddEntry(hFractionProtons, "p", "p");
  if (plotIdentifiedSpectra)
    legend->AddEntry(hFractionProtonsMC, "p", "p");
  legend->AddEntry(hFractionElectrons, "e", "p");
  if (plotIdentifiedSpectra)
    legend->AddEntry(hFractionElectronsMC, "e", "p");
  legend->Draw();
  
  
  
  
  
  
  TH1F* hFractionComparisonPionsPlusMuons = new TH1F(*hFractionPionsPlusMuons);
  hFractionComparisonPionsPlusMuons->SetName("hFractionComparisonPionsPlusMuons");
  hFractionComparisonPionsPlusMuons->GetYaxis()->SetTitle("Ratio to MC truth");
  
  TH1F* hFractionComparisonElectrons = new TH1F(*hFractionElectrons);
  hFractionComparisonElectrons->SetName("hFractionComparisonElectrons");
  hFractionComparisonElectrons->GetYaxis()->SetTitle("Ratio to MC truth");
  
  TH1F* hFractionComparisonKaons = new TH1F(*hFractionKaons);
  hFractionComparisonKaons->SetName("hFractionComparisonKaons");
  hFractionComparisonKaons->GetYaxis()->SetTitle("Ratio to MC truth");
  
  TH1F* hFractionComparisonProtons = new TH1F(*hFractionProtons);
  hFractionComparisonProtons->SetName("hFractionComparisonProtons");
  hFractionComparisonProtons->GetYaxis()->SetTitle("Ratio to MC truth");
  
  hFractionComparisonPionsPlusMuons->Divide(hFractionPionsMC);
  hFractionComparisonElectrons->Divide(hFractionElectronsMC);
  hFractionComparisonKaons->Divide(hFractionKaonsMC);
  hFractionComparisonProtons->Divide(hFractionProtonsMC);
  
  
  TCanvas* cFractionComparisons = new TCanvas("cFractionComparisons", "Particle fraction comparisons",100,10,1200,800);
  cFractionComparisons->SetGridx(1);
  cFractionComparisons->SetGridy(1);
  hFractionComparisonPionsPlusMuons->GetYaxis()->SetRangeUser(0.75, 1.25);
  hFractionComparisonPionsPlusMuons->GetXaxis()->SetMoreLogLabels(kTRUE);
  hFractionComparisonPionsPlusMuons->GetXaxis()->SetNoExponent(kTRUE);
  hFractionComparisonPionsPlusMuons->Draw("e p");
  
  hFractionComparisonElectrons->Draw("e p same");
  hFractionComparisonKaons->Draw("e p same");
  hFractionComparisonProtons->Draw("e p same");
  
  TLegend* legend = new TLegend(0.622126, 0.605932, 0.862069, 0.855932);    
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetNColumns(2);
  legend->AddEntry(hFractionComparisonPionsPlusMuons, "#pi + #mu", "p");
  legend->AddEntry(hFractionComparisonKaons, "K", "p");
  legend->AddEntry(hFractionComparisonProtons, "p", "p");
  legend->AddEntry(hFractionComparisonElectrons, "e", "p");
  legend->Draw();
}



// Xianguo's style
TCanvas* getCanvasLikeXianguo() {
  TCanvas *cc = new TCanvas("cc", "",1,47,600,500);
  cc->SetHighLightColor(2);
  cc->Range(0,0,1,1);
  cc->SetFillColor(0);
  cc->SetBorderMode(0);
  cc->SetBorderSize(2);
  cc->SetFrameBorderMode(0);
  
  return cc;
}

TPad* getPadLikeXianguo(Bool_t upperPanel) {
  TPad* c = 0x0;
  
  if (upperPanel) {
    c = new TPad("c0", "",0,0.4,1,1);
    c->Range(-0.1820482,-0.04943005,1.046867,0.8935751);
    c->SetFillColor(0);
    c->SetBorderMode(0);
    c->SetBorderSize(2);
    c->SetTickx(1);
    c->SetTicky(1);
    c->SetGrid(0,0);
    c->SetLeftMargin(0.14);
    c->SetRightMargin(0.03);
    c->SetTopMargin(0.025);
    c->SetBottomMargin(0.01);
    c->SetFrameBorderMode(0);
    c->SetFrameBorderMode(0);
  }
  else {
    c = new TPad("c1", "",0,0,1,0.4);
    c->Range(-0.1820482,0.14,1.046867,1.473333);
    c->SetFillColor(0);
    c->SetBorderMode(0);
    c->SetBorderSize(2);
    c->SetTickx(1);
    c->SetTicky(1);
    c->SetGrid(0,0);
    c->SetLeftMargin(0.14);
    c->SetRightMargin(0.03);
    c->SetTopMargin(0.01);
    c->SetBottomMargin(0.3);
    c->SetFrameBorderMode(0);
    c->SetFrameBorderMode(0);
  }
  
  return c;
}

void setLegendStyleLikeXianguo(TLegend* leg) {
  leg->SetBorderSize(-1);
  leg->SetTextSize(0.056);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(19);
  leg->SetX1NDC(leg->GetX1NDC() + 0.005);
  leg->SetX2NDC(leg->GetX2NDC() + 0.005);
  leg->SetY1NDC(leg->GetY1NDC() - 0.055);
  leg->SetY2NDC(leg->GetY2NDC() - 0.02);
}

void setLatexStyleLikeXianguo(TLatex* lt) {
  const Double_t textSize = 0.056*2./3.;
  lt->SetTextSize(textSize);
  lt->SetX(lt->GetX() - 0.03);
  lt->SetY(lt->GetY() + 0.03);
}

void setSizesAndOffsetsForHistoLikeXianguo(TH1* h, Bool_t upperPanel) {
  setSizesAndOffsetsLikeXianguo(h->GetXaxis(), h->GetYaxis(), upperPanel);
}

void setSizesAndOffsetsForGraphLikeXianguo(TGraph* h, Bool_t upperPanel) {
  setSizesAndOffsetsLikeXianguo(h->GetXaxis(), h->GetYaxis(), upperPanel);
}
  
void setSizesAndOffsetsLikeXianguo(TAxis* xAxis, TAxis* yAxis, Bool_t upperPanel) {
  if (upperPanel) {
    xAxis->SetLabelFont(42);
    xAxis->SetLabelOffset(0.01);
    xAxis->SetLabelSize(0.07);
    xAxis->SetTitleSize(0.06);
    xAxis->SetTickLength(0.02);
    xAxis->SetTitleOffset(1.35);
    xAxis->SetTitleFont(42);
    yAxis->SetLabelFont(42);
    yAxis->SetLabelOffset(0.01);
    yAxis->SetLabelSize(0.07);
    yAxis->SetTitleSize(0.06);
    yAxis->SetTickLength(0.02);
    yAxis->SetTitleOffset(1.1);
    yAxis->SetTitleFont(42);
  }
  else {
   xAxis->SetMoreLogLabels();
   xAxis->SetLabelFont(42);
   xAxis->SetLabelOffset(0.01);
   xAxis->SetLabelSize(0.1);
   xAxis->SetTitleSize(0.1);
   xAxis->SetTickLength(0.02);
   xAxis->SetTitleOffset(1.35);
   xAxis->SetTitleFont(42);
   yAxis->SetNdivisions(505);
   yAxis->SetLabelFont(42);
   yAxis->SetLabelOffset(0.01);
   yAxis->SetLabelSize(0.1);
   yAxis->SetTitleSize(0.1);
   yAxis->SetTickLength(0.02);
   yAxis->SetTitleOffset(0.6);
   yAxis->SetTitleFont(42);
  }
}


TCanvas* convertCanvasToXianguoStyle(TCanvas* cOld, Bool_t setNewRanges) {
  TCanvas* cNew = getCanvasLikeXianguo();
  TPad* c0 = getPadLikeXianguo(kTRUE);
  TPad* c1 = getPadLikeXianguo(kFALSE);
  
  cNew->cd();
  c0->Draw();
  c1->Draw();
  
  TPad* c0old = ((TPad*)(cOld->GetListOfPrimitives()->At(0)));
  TPad* c1old = ((TPad*)(cOld->GetListOfPrimitives()->At(1)));
  
  for (Int_t i = 0; i < c0old->GetListOfPrimitives()->GetEntries(); i++) {
    if (c0old->GetListOfPrimitives()->At(i)->InheritsFrom("TH1")) {
      TH1D* hd = dynamic_cast<TH1D*>(c0old->GetListOfPrimitives()->At(i));
      setSizesAndOffsetsForHistoLikeXianguo(hd, kTRUE);
      c0->cd();
      if (setNewRanges) {
        hd->GetYaxis()->SetRangeUser(-0.02, 0.94);
        hd->GetXaxis()->SetNdivisions(406);
      }
      TH1D* hDummy = new TH1D(*hd);
      hDummy->SetDirectory(0);
      hDummy->SetName("hDummy");
      for (Int_t i = 1; i <= hDummy->GetNbinsX(); i++) {
        hDummy->SetBinContent(i, -999);
        hDummy->SetBinError(i, 1);
      }
      hDummy->Draw("E1");
      break;//Only draw first hist
    }
  }
  
  for (Int_t i = 0; i < c0old->GetListOfPrimitives()->GetEntries(); i++) {
    if (c0old->GetListOfPrimitives()->At(i)->InheritsFrom("TGraph")) {
      TGraph* g = dynamic_cast<TGraph*>(c0old->GetListOfPrimitives()->At(i));
      setSizesAndOffsetsForGraphLikeXianguo(g, kTRUE);
      c0->cd();
      TString grName(g->GetName());
      TString drawCmd = "E2 psame";
      if (grName.Contains("grstat"))
        drawCmd = "E1 psame";
      g->Draw(drawCmd.Data());
    }
  }
  
  for (Int_t i = 0; i < c0old->GetListOfPrimitives()->GetEntries(); i++) {
    if (c0old->GetListOfPrimitives()->At(i)->InheritsFrom("TLegend")) {
      TLegend* leg = dynamic_cast<TLegend*>(c0old->GetListOfPrimitives()->At(i));
      c0->cd();
      setLegendStyleLikeXianguo(leg);
      leg->Draw("");
      continue;
    }
    
    if (c0old->GetListOfPrimitives()->At(i)->InheritsFrom("TH1")) {
      TH1* h = dynamic_cast<TH1*>(c0old->GetListOfPrimitives()->At(i));
      setSizesAndOffsetsForHistoLikeXianguo(h, kTRUE);
      c0old->cd(); //Must cd to pad in which object was drawn in order to get correct draw option
      const Option_t* opt = h->GetDrawOption();
      c0->cd();
      h->SetDrawOption(opt);
      h->Draw("E1 same");
    }
  }
  
  for (Int_t i = 0; i < cOld->GetListOfPrimitives()->GetEntries(); i++) {
    if (cOld->GetListOfPrimitives()->At(i)->InheritsFrom("TLatex")) {
      TLatex* lt = dynamic_cast<TLatex*>(cOld->GetListOfPrimitives()->At(i));
      setLatexStyleLikeXianguo(lt);
      cNew->cd();
      lt->Draw("");
      break;
    }
  }
  
  c0->Update();    // Update in order to have access to the title in the following
  TPaveText* paveTextTitle = (TPaveText*)c0->FindObject("title");
  if (paveTextTitle) 
    paveTextTitle->Clear();
  c0->Modified();
  
  
  for (Int_t i = 0; i < c1old->GetListOfPrimitives()->GetEntries(); i++) {
    if (c1old->GetListOfPrimitives()->At(i)->InheritsFrom("TH1")) {
      TH1D* hd = dynamic_cast<TH1D*>(c1old->GetListOfPrimitives()->At(i));
      setSizesAndOffsetsForHistoLikeXianguo(hd, kFALSE);
      c1->cd();
      if (setNewRanges) {
        hd->GetYaxis()->SetRangeUser(0.54, 1.46);
        hd->GetXaxis()->SetNdivisions(406);
      }
      TH1D* hDummy = new TH1D(*hd);
      hDummy->SetDirectory(0);
      hDummy->SetName("hDummy");
      for (Int_t i = 1; i <= hDummy->GetNbinsX(); i++) {
        hDummy->SetBinContent(i, -999);
        hDummy->SetBinError(i, 1);
      }
      hDummy->Draw("E1");
      break;//Only draw first hist
    }
  }
  
  
  
  for (Int_t i = 0; i < c1old->GetListOfPrimitives()->GetEntries(); i++) {
    if (c1old->GetListOfPrimitives()->At(i)->InheritsFrom("TGraph")) {
      TGraph* g = dynamic_cast<TGraph*>(c1old->GetListOfPrimitives()->At(i));
      setSizesAndOffsetsForGraphLikeXianguo(g, kTRUE);
      c1->cd();
      TString grName(g->GetName());
      TString drawCmd = "E2 psame";
      if (grName.Contains("grstat"))
        drawCmd = "E1 psame";
      g->Draw(drawCmd.Data());
    }
  }
  
  for (Int_t i = 0; i < c1old->GetListOfPrimitives()->GetEntries(); i++) {
    if (!c1old->GetListOfPrimitives()->At(i)->InheritsFrom("TH1"))
      continue;
    
    TH1* h = dynamic_cast<TH1*>(c1old->GetListOfPrimitives()->At(i));
    setSizesAndOffsetsForHistoLikeXianguo(h, kFALSE);
    c1old->cd(); //Must cd to pad in which object was drawn in order to get correct draw option
    const Option_t* opt = h->GetDrawOption();
    c1->cd();
    h->SetDrawOption(opt);
    h->Draw("E1 same");
  }
  
  c1->Update();    // Update in order to have access to the title in the following
  paveTextTitle = (TPaveText*)c1->FindObject("title");
  if (paveTextTitle) 
    paveTextTitle->Clear();
  c1->Modified();
  
  cNew->Modified();
  cNew->Update();
  
  
  return cNew;
}