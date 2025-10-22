//-----------
//
// Análise de dados para determinar o espectro de energia de um feixe de muões
//
// Hardcoded:
// - número de muões gerados
// - binning dos histogramas
//
//-----------

// número de muões gerados
const int Nevents = 1000000; // quanto maior o número de eventos, mais a curva experimental se aproxima da teórica

// constantes
const double Pi = TMath::Pi();
const double rad_to_deg = 180./Pi;
const double m_to_cm = 100.;

// parâmetros físicos
const double dEdX = 2.;                 // Stopping Power em MeV/(g/cm^2) Constante!!!!!!!!!!
const double densityWater = 1.;         // em g/cm^3
const double densityEarth = 2.38;       // em g/cm^3

// parâmetros geométricos
const double tankHeight = 1.2;          // em m
const double tankRadius = 1.8;          // em m
const double z_AMIGA = -2.3;           // em m

// cálculo da energia perdida pelo muão no tanque, em MeV
double energyLoss(double theta_mu_deg, double phi_mu_deg, double x_out, double y_out, double z_out);  // energia perdida pelo muão no tanque, em MeV - depende de ângulo de incidência do muão (theta e phi), de x_out, y_out e z_out que são as coordenadas do ponto de saída do tanque
// cálculo da energia perdida pelo muão no solo, em MeV
//double energyLoss_2(double theta_mu_deg, double phi_mu_deg, double x_out, double y_out, double z_out, double E_loss_2);

// função principal
void espectroMuoes()
{
    // ficheiro com informação dos muões atmosféricos em Auger
        TFile *fIn = new TFile("Desktop/showerMuons.root", "read");
        TTree *tIn = (TTree*)fIn->Get("showerMuons");
        int entries = tIn->GetEntries();
        if(entries<Nevents)
            return;


// variáveis
    double r, phi, x, y, z,x_AMIGA, y_AMIGA, k;
    double E_mu, E_loss, E_loss_2, E_loss_tot, hE_loss_tot;
    double theta_mu_deg, phi_mu_deg;
    double theta_mu_rad, phi_mu_rad;
    double  dist_solo;
    
// objectos a desenhar
    // energia do muão, em MeV
    double E_mu_min = 0.; double E_mu_max = 10000.; int NbinsE = 100;
    TH1D *hE_mu = new TH1D("hE_mu", "hE_mu", NbinsE, E_mu_min, E_mu_max); // gráfico 1
    TH1D *hE_mu_cumulative = new TH1D("hE_mu_cumulative", "hE_mu_cumulative", NbinsE, E_mu_min, E_mu_max); // gráfico 2
    TH1D *hE_mu_cumulative_AMIGA = new TH1D("hE_mu_cumulative_AMIGA", "hE_mu_cumulative_AMIGA", NbinsE, E_mu_min, E_mu_max); // gráfico 3
    
    // ângulo zenital do muão, em graus
    double theta_mu_min_deg = 0.; double theta_mu_max_deg = Pi/2.*rad_to_deg; int NbinsTheta = (int)theta_mu_max_deg;
    TH1D *hTheta = new TH1D("hTheta", "hTheta", NbinsTheta, theta_mu_min_deg, theta_mu_max_deg); // ângulo que de incidência, gráfico 4
    
    // perfis à superfície
    double Rmin = 0.; double Rmax = tankRadius; int NbinsPos = 20;
    TH2D *hAll = new TH2D("hAll", "hAll", NbinsPos, -Rmax, Rmax, NbinsPos, -Rmax, Rmax);// número de muões gerados
    TH2D *hTransmitted = new TH2D("hTransmitted", "hTransmitted", NbinsPos, -Rmax, Rmax, NbinsPos, -Rmax, Rmax);    // número de muões detectados
    TH2D *hFraction = new TH2D("hFraction", "hFraction", NbinsPos, -Rmax, Rmax, NbinsPos, -Rmax, Rmax);             // fracção de muões detectados
    TProfile2D *hE_loss = new TProfile2D("hEloss", "hEloss", NbinsPos, -Rmax, Rmax, NbinsPos, -Rmax, Rmax);         // valor médio da energia perdida pelos muões
    
    // perfil do plano z=-2.3
    double n_min = -10; double n_max = 10; int Ndiv_plane = 64;
    TH2D *hTransmitted_plane = new TH2D("hTransmitted_plane", "hTransmitted_plane", Ndiv_plane, -n_max, n_max, Ndiv_plane, -n_max, n_max);    // número de muões detectados
    TH2D *hFraction_plane = new TH2D("hFraction_plane", "hFraction_plane", Ndiv_plane, -n_max, n_max, Ndiv_plane, -n_max, n_max);             // fracção de muões detectados
    TProfile2D *hE_loss_plane = new TProfile2D("hE_loss_plane", "hE_loss_plane", Ndiv_plane, -n_max, n_max, Ndiv_plane, -n_max, n_max);         // valor médio da energia perdida pelos muões
    
    // perfil do pad AMIGA
    double div_min = -10; double div_max = 10; int NdivAMIGA = 64;
    TH2D *hTransmittedAMIGA = new TH2D("hTransmittedAMIGA", "hTransmittedAMIGA", NdivAMIGA, -div_max, div_max, NdivAMIGA, -div_max, div_max);    // número de muões detectados
    TH2D *hFractionAMIGA = new TH2D("hFractionAMIGA", "hFractionAMIGA", NdivAMIGA, -div_max, div_max, NdivAMIGA, -div_max, div_max);             // fracção de muões detectados
    TProfile2D *hE_lossAMIGA = new TProfile2D("hE_lossAMIGA", "hE_lossAMIGA", NdivAMIGA, -div_max, div_max, NdivAMIGA, -div_max, div_max);         // valor médio da energia perdida pelos muões
    
    // perfil de amiga em função do perfil à superfície
    TProfile2D *hEloss_AM2 = new TProfile2D("hEloss_AM2", "hEloss_AM2", NbinsPos, -Rmax, Rmax, NbinsPos, -Rmax, Rmax);
    TH2D *hTransmitted_AM2 = new TH2D("hTransmitted_AM2", "hTransmitted_AM2", NbinsPos, -Rmax, Rmax, NbinsPos, -Rmax, Rmax);
    TH2D *hFraction_AM2 = new TH2D("hFraction_AM2", "hFraction_AM2", NbinsPos, -Rmax, Rmax, NbinsPos, -Rmax, Rmax);
    TH2D *hAll_AM2 = new TH2D("hAll_AM2", "hAll_AM2", NbinsPos, -Rmax, Rmax, NbinsPos, -Rmax, Rmax);
    
    // perfil dos RPCs à superfície
    TH2D *hTransmitted_rpc = new TH2D("hTransmitted_rpc", "hTransmitted_rpc", NbinsPos, -Rmax, Rmax, NbinsPos, -Rmax, Rmax);
    TH2D *hFraction_rpc = new TH2D("hFraction_rpc", "hFraction_rpc", NbinsPos, -Rmax, Rmax, NbinsPos, -Rmax, Rmax);
    TProfile2D *hEloss_rpc = new TProfile2D("hEloss_rpc", "hEloss_rpc", NbinsPos, -Rmax, Rmax, NbinsPos, -Rmax, Rmax);
    
// distribuições de probabilidade
    // posição radial do muão à superfície
    TF1 *fdistR = new TF1("fdistR", "x", Rmin, Rmax);

    // energia do muão
    // TF1 *fdistE = new TF1("fdistE", "1", E_mu_min, E_mu_max);   // distribuição uniforme
    
    // TF1 *fdistE = new TF1("fdistE", "x", E_mu_min, E_mu_max);   // distribuição linear
    
    //TF1 *fdistE = new TF1("fdistE", "x^3", E_mu_min, E_mu_max);   // distribuição exponencial
    
    
    
    // ângulo zenital do muão
    TF1 *fdistTheta = new TF1("fdistTheta", "TMath::Sin(x)*TMath::Cos(x)*TMath::Cos(x)", theta_mu_min_deg/rad_to_deg, theta_mu_max_deg/rad_to_deg);
    
// geração de eventos
    for(int i=0; i<Nevents; i++)
    {
        // gerar a posição do muão
        r = fdistR->GetRandom();
        phi = gRandom->Uniform(0., 2.*Pi);
        x = r * TMath::Cos(phi);
        y = r * TMath::Sin(phi);
        z = 0.;
        hAll->Fill(x, y);   // preenchimento do perfil respectivo
        
        
        // gerar a energia do muão
//      E_mu = fdistE->GetRandom(); // a energia do muão depende de uma distribuição
                tIn->GetEntry(i);
                E_mu = tIn->GetLeaf("MuEnergyMeV")->GetValue(0);

        hE_mu->Fill(E_mu);
        
        // gerar a direcção do muão
//        theta_mu_deg = fdistTheta->GetRandom()*rad_to_deg;   // direção do muão
//        phi_mu_deg = gRandom->Uniform(0., 2.*Pi)*rad_to_deg;  // não há direções preferenciais ( apesar do campo magnético)
        theta_mu_deg = tIn->GetLeaf("MuThetaDeg")->GetValue(0);
        phi_mu_deg = tIn->GetLeaf("MuPhiDeg")->GetValue(0);
        hTheta->Fill(theta_mu_deg);
        
        // energia perdida pelo muão dentro do tanque
        E_loss = energyLoss(theta_mu_deg, phi_mu_deg, x, y, z); // ponto de saída e a trajetória (cálculo da distância do muão no tanque), Stopping Power 2 MeV/(g/cm^2)
        hE_loss->Fill(x, y, E_loss);
        
        
        
        // 1 detetor: x_AMIGA em [minXa, maxXa] e y_AMIGA em [minYa, maxYa]
            if(x >= -1.4 && x <= -0.2 && y >= 0.2 && y <= 1.7)
            {
                hEloss_rpc->Fill(x, y, E_loss);
            }

        // 2 detetor: x_AMIGA em [minXa, maxXa] e y_AMIGA em [minYa, maxYa]
            else if(x >= 0.2 && x <= 1.4 &&  y >= -1.7 && y <= -0.2)
            {
                hEloss_rpc->Fill(x, y, E_loss);
            }

            
        // 3 detetor: x_AMIGA em [minXa, maxXa] e y_AMIGA em [minYa, maxYa]
            else if(x >= -1.4 && x <= -0.2 &&  y >= -1.7 && y <= -0.2)
            {
                hEloss_rpc->Fill(x, y, E_loss);
            }

        // 4 detetor: x_AMIGA em [minXa, maxXa] e y_AMIGA em [minYa, maxYa]
            else if(x >= 0.2 && x <= 1.4 && y >= 0.2 && y <= 1.7)
            {
                hEloss_rpc->Fill(x, y, E_loss);
            }

        
        
        
//Interseção como plano z=-2.3
        
        k = z_AMIGA / (TMath::Cos(theta_mu_deg/rad_to_deg));
        x_AMIGA = x + k * (TMath::Cos(phi_mu_deg/rad_to_deg)*TMath::Sin(theta_mu_deg/rad_to_deg));
        y_AMIGA = y + k * (TMath::Sin(phi_mu_deg/rad_to_deg)*TMath::Sin(theta_mu_deg/rad_to_deg));
        
//        cout << x_AMIGA <<" "<< y_AMIGA <<" "<< k << " " << z_AMIGA << endl;
    
        
        // comprimento do trajecto do muão no solo até chegar ao detetor AMIGA
            dist_solo = TMath::Sqrt( (x_AMIGA-x)*(x_AMIGA-x) + (y_AMIGA-y)*(y_AMIGA-y) + (z_AMIGA-z)*(z_AMIGA-z));
                
        // energia perdida pelo muão na terra depois de atravessar o tanque
            E_loss_2 = dist_solo*m_to_cm * densityEarth * dEdX;
        
        
        
        E_loss_tot = E_loss+E_loss_2;
        
        // caso o muão não seja absorvido pelo tanque, o mesmo vai ser transmitido (passa do tanque para o solo e, caso tenha uma determinada energia, chega até ao pad "AMIGA").
         hE_loss_plane->Fill(x_AMIGA, y_AMIGA, E_loss_tot); // energia média  perdida numa célula
        
      
        
        
        // condição A: x_AMIGA em [minXa, maxXa] e y_AMIGA em [minYa, maxYa]
            if(x_AMIGA >= -7 && x_AMIGA <= -6 && y_AMIGA >= -8 && y_AMIGA <= 2)
            {
                hE_lossAMIGA->Fill(x_AMIGA, y_AMIGA, E_loss_tot);
                hEloss_AM2->Fill(x, y,E_loss_tot);
                hAll_AM2->Fill(x, y);
                if(E_loss_tot<E_mu){
                    hTransmitted_AM2->Fill(x,y);
                }
               
            }

        // condição B: x_AMIGA em [minXa, maxXa] e y_AMIGA em [minYa, maxYa]
            else if(x_AMIGA >= -6 && x_AMIGA <= 4 && y_AMIGA >= -8 && y_AMIGA <= -7)
            {
                hE_lossAMIGA->Fill(x_AMIGA, y_AMIGA, E_loss_tot);
                hEloss_AM2->Fill(x, y,E_loss_tot);
                hAll_AM2->Fill(x, y);
                if(E_loss_tot<E_mu){
                    hTransmitted_AM2->Fill(x,y);
                }
            }

            
        // condição C: x_AMIGA em [minXa, maxXa] e y_AMIGA em [minYa, maxYa]
            else if(x_AMIGA >= -5 && x_AMIGA <= 0 && y_AMIGA >= -7 && y_AMIGA <= -6)
            {
                hE_lossAMIGA->Fill(x_AMIGA, y_AMIGA, E_loss_tot);
                hEloss_AM2->Fill(x, y,E_loss_tot);
                hAll_AM2->Fill(x, y);
                if(E_loss_tot<E_mu){
                    hTransmitted_AM2->Fill(x,y);
                }
            }

        // condição D: x_AMIGA em [minXa, maxXa] e y_AMIGA em [minYa, maxYa]
            else if(x_AMIGA >= -6 && x_AMIGA <= -5 && y_AMIGA >= -7 && y_AMIGA <= -2)
            {
                hE_lossAMIGA->Fill(x_AMIGA, y_AMIGA, E_loss_tot);
                hEloss_AM2->Fill(x, y,E_loss_tot);
                hAll_AM2->Fill(x, y);
                if(E_loss_tot<E_mu){
                    hTransmitted_AM2->Fill(x,y);
                }
            }

        // condição E: x_AMIGA em [minXa, maxXa] e y_AMIGA em [minYa, maxYa]
            else if(x_AMIGA >= 5 && x_AMIGA <= 6 && y_AMIGA >= -2 && y_AMIGA <= 8)
            {
                hE_lossAMIGA->Fill(x_AMIGA, y_AMIGA, E_loss_tot);
                hEloss_AM2->Fill(x, y,E_loss_tot);
                hAll_AM2->Fill(x, y);
                if(E_loss_tot<E_mu){
                    hTransmitted_AM2->Fill(x,y);
                }
            }

        // condição F: x_AMIGA em [minXa, maxXa] e y_AMIGA em [minYa, maxYa]
            else if(x_AMIGA >= -5 && x_AMIGA <= 5 && y_AMIGA >= 7 && y_AMIGA <= 8)
            {
                hE_lossAMIGA->Fill(x_AMIGA, y_AMIGA, E_loss_tot);
                hEloss_AM2->Fill(x, y,E_loss_tot);
                hAll_AM2->Fill(x, y);
                if(E_loss_tot<E_mu){
                    hTransmitted_AM2->Fill(x,y);
                }
            }

            
        // condição G: x_AMIGA em [minXa, maxXa] e y_AMIGA em [minYa, maxYa]
            else if(x_AMIGA >= -5 && x_AMIGA <= 5 && y_AMIGA >= 6 && y_AMIGA <= 7)
            {
                hE_lossAMIGA->Fill(x_AMIGA, y_AMIGA, E_loss_tot);
                hEloss_AM2->Fill(x, y,E_loss_tot);
                hAll_AM2->Fill(x, y);
                if(E_loss_tot<E_mu){
                    hTransmitted_AM2->Fill(x,y);
                }
            }

        // condição tank: x_AMIGA em [minXa, maxXa] e y_AMIGA em [minYa, maxYa]
//           else if( (x_AMIGA * x_AMIGA) + (y_AMIGA * y_AMIGA) <= (tankRadius * tankRadius))
//            {
//                hE_lossAMIGA->Fill(x_AMIGA, y_AMIGA, E_loss_tot);
//
//            }

        
        // verificar se o muão foi absorvido no tanque
        if(E_loss<E_mu)
        {
            hTransmitted->Fill(x, y);
            
            if(E_loss_tot<E_mu){
                hTransmittedAMIGA->Fill(x_AMIGA, y_AMIGA);
            }
                
        }
    }
    
   
    
// fracção de muões detectados no fundo do tanque
    hFraction->Add(hTransmitted);
    hFraction->Divide(hAll);
    
    
// fracção de muões detectados no fundo do tanque
    hFraction_AM2->Add(hTransmitted_AM2);
    hFraction_AM2->Divide(hAll_AM2);
    
   
    
// espectro de energia dos muões acumulado (linha azul)
    for(int i=1; i<=NbinsE; i++)
        hE_mu_cumulative->SetBinContent(i, hE_mu->Integral(0, i));
    hE_mu_cumulative->Scale(1./(double)Nevents);
    

// fracção de muões absorvidos em função do valor médio da energia perdida pelos muões (pontos pretos)
    TGraph *gr = new TGraph();
    int Ncells = hE_loss->GetNcells();
    for(int i=0; i<Ncells; i++)
        gr->SetPoint(i, hE_loss->GetBinContent(i), 1.-hFraction->GetBinContent(i)); // fração  de eventos absorvidos no tanque
    for(int i=0; i<Ncells; i++)
        gr->SetPoint(i+Ncells, hEloss_AM2->GetBinContent(i), 1.-hFraction_AM2->GetBinContent(i));

// formatos
    //gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    
// gráficos
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 800);
    c1->Divide(4, 4);
    
    c1->cd(1); hE_mu->Draw(); hE_mu->GetXaxis()->SetTitle("E_{#mu} (MeV)"); hE_mu->GetYaxis()->SetTitle("Nevents");
    c1->cd(2); hE_mu_cumulative->Draw(); hE_mu_cumulative->GetXaxis()->SetTitle("E_{#mu} (MeV)"); hE_mu_cumulative->GetYaxis()->SetTitle("#int_{0}^{E} p(x)dx ");
    c1->cd(3); hE_mu_cumulative_AMIGA->Draw(); hE_mu_cumulative_AMIGA->GetXaxis()->SetTitle("E_{#mu} (MeV)"); hE_mu_cumulative_AMIGA->GetYaxis()->SetTitle("#int_{0}^{E} p(x)dx ");
    c1->cd(4); hTheta->Draw(); hTheta->GetXaxis()->SetTitle("#theta_{#mu} (#circ)"); hTheta->GetYaxis()->SetTitle("Nevents");
    c1->cd(5); hE_loss->Draw("zcol"); hE_loss->GetXaxis()->SetTitle("x (m)"); hE_loss->GetYaxis()->SetTitle("y (m)"); hE_loss->SetTitle("Energy loss tank base (MeV)");
    c1->cd(6); hEloss_rpc->Draw("zcol"); hEloss_rpc->GetXaxis()->SetTitle("x (m)"); hEloss_rpc->GetYaxis()->SetTitle("y (m)"); hEloss_rpc->SetTitle("RPCs Energy loss (MeV)");
    c1->cd(7); hFraction->Draw("zcol"); hFraction->GetXaxis()->SetTitle("x (m)"); hFraction->GetYaxis()->SetTitle("y (m)"); hFraction->SetTitle("Transmitted fraction");
    c1->cd(8); hE_loss_plane->Draw("zcol"); hE_loss_plane->GetXaxis()->SetTitle("x (m)"); hE_loss_plane->GetYaxis()->SetTitle("y (m)"); hE_loss_plane->SetTitle("Energy loss z= -2.3 (MeV)");
    c1->cd(9); hE_lossAMIGA->Draw("zcol"); hE_lossAMIGA->GetXaxis()->SetTitle("x (m)"); hE_lossAMIGA->GetYaxis()->SetTitle("y (m)"); hE_lossAMIGA->SetTitle("AMIGA's Energy loss");
    c1->cd(10); hEloss_AM2->Draw("zcol"); hEloss_AM2->GetXaxis()->SetTitle("x (m)"); hEloss_AM2->GetYaxis()->SetTitle("y (m)"); hEloss_AM2->SetTitle("");
    c1->cd(11); hFraction_AM2->Draw("zcol"); hFraction_AM2->GetXaxis()->SetTitle("x (m)"); hFraction_AM2->GetYaxis()->SetTitle("y (m)"); hFraction_AM2->SetTitle("");
    
// ajuste de uma função aos "pontos experimentais" tank (linha vermelha)
    c1->cd(2);
    TF1 *ferf = new TF1("ferf", "", E_mu_min, E_mu_max);
    ferf->SetParameter(0, 1);
    ferf->SetParameter(1, 0);
    ferf->FixParameter(2, 0);
//    gr->Fit(ferf,"","",250.,300.);
    gr->Draw("psame");
    ferf->Draw("same");
    

}



// cálculo da energia perdida pelo muão no tanque, em MeV
double energyLoss(double theta_mu_deg, double phi_mu_deg, double x_out, double y_out, double z_out)
{
    double x_in, y_in, z_in;
    double k0, R, a, b, c, arg, k1, k2;
    double dist, E_loss;
    
    double theta_mu_rad = theta_mu_deg/rad_to_deg;
    double sinTheta_mu = TMath::Sin(theta_mu_rad);
    double cosTheta_mu = TMath::Cos(theta_mu_rad);
    
    double phi_mu_rad = phi_mu_deg/rad_to_deg;
    double sinPhi_mu = TMath::Sin(phi_mu_rad);
    double cosPhi_mu = TMath::Cos(phi_mu_rad);
    
// ponto de entrada do muão no tanque
    
    // posição radial do muão para z=tankHeight
    k0 = (tankHeight-z_out) / cosTheta_mu;
    
    R = x_out*x_out + y_out*y_out;
    R += 2. * k0 * sinTheta_mu * (cosPhi_mu*x_out + sinPhi_mu*y_out);
    R += k0*k0 * sinTheta_mu*sinTheta_mu;
    R = TMath::Sqrt(R);
    
    // o muão entra no tanque pelo topo
    if(R<=tankRadius)
    {
        x_in = x_out + k0 * sinTheta_mu * cosPhi_mu;
        y_in = y_out + k0 * sinTheta_mu * sinPhi_mu;
        z_in = z_out + tankHeight;
    }
    
    // o muão entra no tanque pelo lado
    else
    {
        a = sinTheta_mu*sinTheta_mu;
        b = 2. * sinTheta_mu * (cosPhi_mu*x_out + sinPhi_mu*y_out);
        c = x_out*x_out + y_out*y_out - tankRadius*tankRadius;
        
        arg = b*b - 4.*a*c;
        k1 = -(b-TMath::Sqrt(arg)) / (2.*a);
        k2 = -(b+TMath::Sqrt(arg)) / (2.*a);
        
        x_in = x_out + k1 * sinTheta_mu * cosPhi_mu;
        y_in = y_out + k1 * sinTheta_mu * sinPhi_mu;
        z_in = z_out + k1 * cosTheta_mu;
        
        if(z_in<z_out || z_in>z_out+tankHeight) // condição: pelo menos 1 verdadeira
        {
            x_in = x_out + k2 * sinTheta_mu * cosPhi_mu;
            y_in = y_out + k2 * sinTheta_mu * sinPhi_mu;
            z_in = z_out + k2 * cosTheta_mu;
        }
    }
    
// comprimento do trajecto do muão dentro do tanque
    dist = TMath::Sqrt( (x_out-x_in)*(x_out-x_in) + (y_out-y_in)*(y_out-y_in) + (z_out-z_in)*(z_out-z_in) );
    
// energia perdida pelo muão dentro do tanque
    E_loss = dist*m_to_cm * densityWater * dEdX;
    
    return E_loss;
}
    
