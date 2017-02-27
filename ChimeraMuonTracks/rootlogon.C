{
    //gSystem->Setenv("MCTABLES","../tables");

    //gSystem->Load("libRooFit.so");

    //gSystem->Load("/usr/local/lib/libcfitsio.so");
    //gSystem->Load("/usr/lib/libcfitsio.so");

    // Load RooFit
    //gSystem->Load("libRooFit.so");

    // Load MaxCam libraries
    //  gSystem->Load("/usr/local/lib/MaxCam.so"); // production copy
    //gSystem->Load("../MaxCam_linux.so"); // local copy
    //gSystem->Load("../MaxCam_macosx.so"); // local copy


    gStyle->SetPalette(1);
    //TGaxis::SetMaxDigits(3);

    gStyle->SetFrameBorderMode(0);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetPadColor(0);
    gStyle->SetCanvasColor(0);
    gStyle->SetTitleColor(1);
    gStyle->SetTitleFont(42);
    gStyle->SetStatColor(1);
    //gStyle->SetFillColor(0);
    gStyle->SetTitleColor(1);

    // set the paper & margin sizes
    //gStyle->SetPaperSize(20,26);
    gStyle->SetPadTopMargin(0.13);
    gStyle->SetPadRightMargin(0.13);
    gStyle->SetPadBottomMargin(0.13);
    gStyle->SetPadLeftMargin(0.13);

    gStyle->SetNdivisions(1010,"x");
    gStyle->SetNdivisions(1010,"y");

    // use large Times-Roman fonts
    gStyle->SetTextFont(42);
    gStyle->SetTextSize(0.08);

    gStyle->SetLabelFont(42,"x");
    gStyle->SetLabelFont(42,"y");
    gStyle->SetLabelFont(42,"z");
    gStyle->SetTitleFont(42,"x");
    gStyle->SetTitleFont(42,"y");
    gStyle->SetTitleFont(42,"z");

    gStyle->SetTitleOffset(1.00, "x");
    gStyle->SetTitleOffset(1.15, "y");
    gStyle->SetTitleOffset(1.00, "z");

    gStyle->SetLabelSize(0.05,"x");
    gStyle->SetTitleSize(0.05,"x");
    gStyle->SetLabelSize(0.05,"y");
    gStyle->SetTitleSize(0.05,"y");
    gStyle->SetLabelSize(0.05,"z");
    gStyle->SetTitleSize(0.05,"z");

    // stat box
    gStyle->SetStatBorderSize(1);
    gStyle->SetStatX(1-gStyle->GetPadRightMargin());
    gStyle->SetStatY(1-gStyle->GetPadTopMargin());
    gStyle->SetStatW(.200);
    gStyle->SetStatH(.125);
    gStyle->SetStatColor(0);
    gStyle->SetStatStyle(1001);
    gStyle->SetStatFont(42);

    // title
    gStyle->SetTitleX(0.5f);
    gStyle->SetTitleW(0.8f);


    // use bold lines and markers
    gStyle->SetMarkerStyle(1);
    gStyle->SetHistLineWidth(2);
    gStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

    // do not display any of the standard histogram decorations
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    
    // put tick marks on top and RHS of plots
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetGridStyle(3);
    gStyle->SetPadGridX(true);
    gStyle->SetPadGridY(true);
    
    
    gStyle->SetPalette(1);
}

