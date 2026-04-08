Comparación codigo Santiago y mi código.
Mi codigo es el llamado ExtractVariables_Etapa2.C

En este archivo voy a poner como santiago define las cosas para poder comparar sin tener que leer todo el archivo. De todas formas el archivo de Santiago se llama FSRGammaGammaHHbbbb.C y esta en el repositorio XCC_HH_bbb.

Solo cambian cosas en ExtractVariables_Etapa2.C porque es donde se definen las funciones. El resto de las cosas se sacan directo del output de delphes.

---

## Forma del evento (Sphericity o Aplanarity)

```cpp
float findEventShape(TLorentzVector j1, TLorentzVector j2, TLorentzVector j3, TLorentzVector j4, string shape)
{
	vector<TVector3> vectorCollection = {
		TVector3(j1.Px(), j1.Py(), j1.Pz()),
		TVector3(j2.Px(), j2.Py(), j2.Pz()),
		TVector3(j3.Px(), j3.Py(), j3.Pz()),
		TVector3(j4.Px(), j4.Py(), j4.Pz())
	};
	TMatrixDSym momentumTensor(3);
	for(std::vector<TVector3>::const_iterator p=vectorCollection.begin(); p!=vectorCollection.end(); p++){
		for(int k=0; k<3; k++){
			for(int m=0; m<=k; m++){
				momentumTensor[k][m] += (*p)[k]*(*p)[m];
			}
		}
	}
	momentumTensor *= 1/(momentumTensor[0][0] + momentumTensor[1][1] + momentumTensor[2][2]);
	TMatrixDSymEigen eigenSystem(momentumTensor);
	TVectorD eigenvalues = eigenSystem.GetEigenValues();
	float eigenvalue1_ = eigenvalues[0];
	float eigenvalue2_ = eigenvalues[1];
	float eigenvalue3_ = eigenvalues[2];
	if(shape == "sphericity") return 1.5*(eigenvalue2_ + eigenvalue3_)/(eigenvalue1_ + eigenvalue2_ + eigenvalue3_);
	else if(shape == "aplanarity") return 1.5*(eigenvalue3_)/(eigenvalue1_ + eigenvalue2_ + eigenvalue3_);
	else return 0;
}
```

---

## Thrust

```cpp
float findThrust(vector<TLorentzVector>& momenta, TLorentzVector& thrustAxis)
{
	float maxThrust = 0.0;
	for(float theta = 0.0; theta < TMath::Pi(); theta += 0.1)
	{
		for(float phi = 0.0; phi < 2.0 * TMath::Pi(); phi += 0.1)
		{
			TLorentzVector axis(TMath::Sin(theta) * TMath::Cos(phi),
			                    TMath::Sin(theta) * TMath::Sin(phi),
			                    TMath::Cos(theta), 0.0);
			float sumP = 0.0;
			float sumP_parallel = 0.0;
			for (auto& p : momenta)
			{
				sumP += p.P();
				sumP_parallel += TMath::Abs(p.Dot(axis) / Magnitude(axis));
			}
			float thrust = sumP_parallel / sumP;
			if (thrust > maxThrust)
			{
				maxThrust = thrust;
				thrustAxis = axis;
			}
		}
	}
	return maxThrust;
}
```

---

## JetPairs HH

Calcula los chiSquared ZZ pero no los usa (codigo de mas). Tambien recalcula los pares dentro del if/else cuando ya los habia calculado antes.

```cpp
void findJetPairs(TLorentzVector& jetPairB1, TLorentzVector& jetPairB2,
                  TLorentzVector jetB1, TLorentzVector jetB2,
                  TLorentzVector jetB3, TLorentzVector jetB4,
                  double& jetPairB1Index1, double& jetPairB1Index2,
                  double& jetPairB2Index1, double& jetPairB2Index2,
                  double& minJetChiS)
{
	double errorJetPairMass=1;
	jetPairB1=jetB1+jetB2;
	jetPairB2=jetB3+jetB4;
	double jetChiS12=pow((jetPairB1.M()-125), 2)/errorJetPairMass + pow((jetPairB2.M()-125), 2)/errorJetPairMass;
	double jetChiS12ZZ=pow((jetPairB1.M()-90), 2)/errorJetPairMass + pow((jetPairB2.M()-90), 2)/errorJetPairMass; // nunca se usa
	jetPairB1=jetB1+jetB3;
	jetPairB2=jetB2+jetB4;
	double jetChiS13=pow((jetPairB1.M()-125), 2)/errorJetPairMass + pow((jetPairB2.M()-125), 2)/errorJetPairMass;
	double jetChiS13ZZ=pow((jetPairB1.M()-90), 2)/errorJetPairMass + pow((jetPairB2.M()-90), 2)/errorJetPairMass; // nunca se usa
	jetPairB1=jetB1+jetB4;
	jetPairB2=jetB2+jetB3;
	double jetChiS14=pow((jetPairB1.M()-125), 2)/errorJetPairMass + pow((jetPairB2.M()-125), 2)/errorJetPairMass;
	minJetChiS=TMath::Min(TMath::Min(jetChiS12, jetChiS13), jetChiS14);
	if(minJetChiS==jetChiS12)
	{
		jetPairB1=jetB1+jetB2; // redundante, ya calculado
		jetPairB2=jetB3+jetB4;
		jetPairB1Index1=0; jetPairB1Index2=1;
		jetPairB2Index1=2; jetPairB2Index2=3;
	}
	else if(minJetChiS==jetChiS13)
	{
		jetPairB1=jetB1+jetB3;
		jetPairB2=jetB2+jetB4;
		jetPairB1Index1=0; jetPairB1Index2=2;
		jetPairB2Index1=1; jetPairB2Index2=3;
	}
	else if(minJetChiS==jetChiS14)
	{
		jetPairB1=jetB1+jetB4;
		jetPairB2=jetB2+jetB3;
		jetPairB1Index1=0; jetPairB1Index2=3;
		jetPairB2Index1=1; jetPairB2Index2=2;
	}
}
```

---

## JetPairs ZZ

Tiene contadores de ventanas en masa Z (probablemente para debugging) y recalcula los pares redundantemente igual que en HH.

```cpp
void findJetPairsZZ(TLorentzVector& jetPairB1, TLorentzVector& jetPairB2,
                    TLorentzVector jetB1, TLorentzVector jetB2,
                    TLorentzVector jetB3, TLorentzVector jetB4,
                    double& jetPairB1Index1, double& jetPairB1Index2,
                    double& jetPairB2Index1, double& jetPairB2Index2,
                    bool& flagZZMass, float& minChiSquaredZZMass,
                    double distanceZMass,
                    int& contEventsZWindow10, int& contEventsZWindow01,
                    int& contEventsZWindow05, int& contEventsZWindow1,
                    int& contEventsZWindow15, int& contEventsZWindow2,
                    int& contEventsZWindow5,
                    float& distanceZ1MinChiSquaredZZMass,
                    float& distanceZ2MinChiSquaredZZMass)
{
	double errorJetPairMass=1;
	jetPairB1=jetB1+jetB2;
	jetPairB2=jetB3+jetB4;
	double jetChiS12=pow((jetPairB1.M()-90), 2)/errorJetPairMass + pow((jetPairB2.M()-90), 2)/errorJetPairMass;
	jetPairB1=jetB1+jetB3;
	jetPairB2=jetB2+jetB4;
	double jetChiS13=pow((jetPairB1.M()-90), 2)/errorJetPairMass + pow((jetPairB2.M()-90), 2)/errorJetPairMass;
	jetPairB1=jetB1+jetB4;
	jetPairB2=jetB2+jetB3;
	double jetChiS14=pow((jetPairB1.M()-90), 2)/errorJetPairMass + pow((jetPairB2.M()-90), 2)/errorJetPairMass;
	minChiSquaredZZMass = TMath::Min(TMath::Min(jetChiS12, jetChiS13), jetChiS14);
	if(minChiSquaredZZMass==jetChiS12)
	{
		jetPairB1=jetB1+jetB2; // redundante
		jetPairB2=jetB3+jetB4;
		jetPairB1Index1=0; jetPairB1Index2=1;
		jetPairB2Index1=2; jetPairB2Index2=3;
	}
	else if(minChiSquaredZZMass==jetChiS13)
	{
		jetPairB1=jetB1+jetB3;
		jetPairB2=jetB2+jetB4;
		jetPairB1Index1=0; jetPairB1Index2=2;
		jetPairB2Index1=1; jetPairB2Index2=3;
	}
	else if(minChiSquaredZZMass==jetChiS14)
	{
		jetPairB1=jetB1+jetB4;
		jetPairB2=jetB2+jetB3;
		jetPairB1Index1=0; jetPairB1Index2=3;
		jetPairB2Index1=1; jetPairB2Index2=2;
	}
	if(abs(jetPairB1.M()-90)<distanceZMass && abs(jetPairB2.M()-90)<distanceZMass) flagZZMass=true;
	if(abs(jetPairB1.M()-90)<0.1  && abs(jetPairB2.M()-90)<0.1)  contEventsZWindow01++;  // contadores para debugging
	if(abs(jetPairB1.M()-90)<0.5  && abs(jetPairB2.M()-90)<0.5)  contEventsZWindow05++;
	if(abs(jetPairB1.M()-90)<1    && abs(jetPairB2.M()-90)<1)    contEventsZWindow1++;
	if(abs(jetPairB1.M()-90)<1.5  && abs(jetPairB2.M()-90)<1.5)  contEventsZWindow15++;
	if(abs(jetPairB1.M()-90)<2    && abs(jetPairB2.M()-90)<2)    contEventsZWindow2++;
	if(abs(jetPairB1.M()-90)<5    && abs(jetPairB2.M()-90)<5)    contEventsZWindow5++;
	if(abs(jetPairB1.M()-90)<10   && abs(jetPairB2.M()-90)<10)   contEventsZWindow10++;
	distanceZ1MinChiSquaredZZMass=abs(jetPairB1.M()-90);
	distanceZ2MinChiSquaredZZMass=abs(jetPairB2.M()-90);
}
```

---

## MinJetM

Guarda los 6 TLorentzVector en variables nombradas antes de calcular la masa. Funcionalmente identico a mi version.

```cpp
float findMinJetM(TLorentzVector jetB1, TLorentzVector jetB2,
                  TLorentzVector jetB3, TLorentzVector jetB4)
{
	TLorentzVector jetPairBB1=jetB1+jetB2;
	TLorentzVector jetPairBB2=jetB3+jetB4;
	TLorentzVector jetPairBB3=jetB1+jetB3;
	TLorentzVector jetPairBB4=jetB2+jetB4;
	TLorentzVector jetPairBB5=jetB1+jetB4;
	TLorentzVector jetPairBB6=jetB2+jetB3;
	float minMass = TMath::Min(static_cast<float>(jetPairBB1.M()), static_cast<float>(jetPairBB2.M()));
	minMass = TMath::Min(minMass, static_cast<float>(jetPairBB3.M()));
	minMass = TMath::Min(minMass, static_cast<float>(jetPairBB4.M()));
	minMass = TMath::Min(minMass, static_cast<float>(jetPairBB5.M()));
	minMass = TMath::Min(minMass, static_cast<float>(jetPairBB6.M()));
	return minMass;
}
```
