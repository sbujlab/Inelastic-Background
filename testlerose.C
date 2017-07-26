#include "leroseres.h"

void testlerose(){
    LoadLeroseRes();

    TH1F *h = new TH1F("h", "h", 400, -0.03, 0.07);

    int i;
    for( i = 0; i < 1000000; i++ ){
	h->Fill(getlerose_ressmear()*2.2);
    }

    h->Draw();
}
