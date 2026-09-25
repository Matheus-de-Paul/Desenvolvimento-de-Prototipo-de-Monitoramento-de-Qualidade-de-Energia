# Desenvolvimento-de-Prototipo-de-Monitoramento-de-Qualidade-de-Energia
/*
 * ================================================================
 *  Protótipo de monitoramento da qualidade de energia - ESP32
 * ================================================================
 *  Tensão:   ZMPT101B -> divisor 10k / 20k -> GPIO34
 *  Corrente: SCT-013 30A/1V -> polarização 1,65V (R4/R5) -> 1k -> GPIO35
 *
 *  Mede: frequência, tensão RMS, corrente RMS, potência ativa,
 *        reativa e aparente, fator de potência, defasagem V-I,
 *        harmônicos (FFT) de tensão e corrente e THD.
 *
 *  Bibliotecas: arduinoFFT (Enrique Condes), versão 2.x
 * ================================================================
 */
#include <arduinoFFT.h>

// ===================== PINOS =====================
const int PINO_TENSAO   = 34;
const int PINO_CORRENTE = 35;

// ===================== AMOSTRAGEM =====================
const int N_AMOSTRAS     = 64;   // amostras por ciclo da rede (potência de 2 para a FFT)
const int N_CICLOS_MEDIA = 30;   // ciclos somados em cada relatório (~0,5 s)
const int N_HARMONICOS   = 10;   // harmônicos calculados (60 Hz até 600 Hz)

// ===================== CALIBRAÇÃO =====================
// FATOR_TENSAO   = tensão do multímetro      / "RMS bruto V"
// FATOR_CORRENTE = corrente do alicate × N_VOLTAS / "RMS bruto I"
const double FATOR_TENSAO   = 0.30;   // V por contagem do ADC
const double FATOR_CORRENTE = 0.024;  // A por contagem do ADC (estimativa inicial)
const int    N_VOLTAS       = 2;      // voltas do fio na pinça (motor: 1)
const double RUIDO_I_RMS    = 0.0;    // corrente lida SEM carga (A), descontada do resultado
const bool   INVERTER_CORRENTE = false;     // true se a potência ativa sair negativa
// Correção de fase do sensor de corrente (graus). Calibre com uma carga resistiva
// (lâmpada incandescente): coloque aqui a "Defasagem V-I" que aparecer com ela.
const double CORRECAO_FASE_GRAUS = 23.6;
const bool   TRIFASICO_EQUILIBRADO = false; // true no motor: mostra totais trifásicos (x3)

// ===================== LIMIARES =====================
const double VRMS_MINIMO = 20.0;  // abaixo disso: considera que não há tensão (V)
const double IRMS_MINIMO = 0.10;  // abaixo disso: considera que não há corrente (A)

// ===================== VARIÁVEIS =====================
uint16_t amostrasV[N_AMOSTRAS];
uint16_t amostrasI[N_AMOSTRAS];
double vReal[N_AMOSTRAS];
double vImag[N_AMOSTRAS];
ArduinoFFT<double> FFT = ArduinoFFT<double>(vReal, vImag, N_AMOSTRAS, 7680.0);

double freqRede    = 60.0;  // frequência usada para sincronizar a amostragem
double mediaVAtual = 2048;  // offset da tensão (usado na medição de frequência)

// Acumuladores do relatório
double somaV2, somaI2, somaP, somaBrutoV2, somaBrutoI2;
double harmV2[N_HARMONICOS + 1];
double harmI2[N_HARMONICOS + 1];
double cruzRe, cruzIm;      // produto cruzado das fundamentais (para a defasagem)
int nCiclos;
int nAtrasos;               // amostras que saíram atrasadas (ADC lento demais)
double duracaoBlocoUs;      // duração real do último ciclo coletado

// ================================================================
void zerarAcumuladores() {
  somaV2 = somaI2 = somaP = somaBrutoV2 = somaBrutoI2 = 0;
  cruzRe = cruzIm = 0;
  for (int h = 0; h <= N_HARMONICOS; h++) { harmV2[h] = 0; harmI2[h] = 0; }
  nCiclos = 0;
  nAtrasos = 0;
}

// Lê o pino duas vezes e faz a média (reduz o ruído do ADC)
uint16_t lerADC(int pino) {
  return (analogRead(pino) + analogRead(pino)) / 2;
}

// Lê uma amostra de tensão e uma de corrente na ordem V, I, I, V.
// Assim as duas médias ficam centradas no mesmo instante (sem atraso entre V e I).
void lerAmostra(uint16_t &v, uint16_t &i) {
  uint32_t sv = analogRead(PINO_TENSAO);
  uint32_t si = analogRead(PINO_CORRENTE);
  si += analogRead(PINO_CORRENTE);
  sv += analogRead(PINO_TENSAO);
  v = sv / 2;
  i = si / 2;
}

// Coleta 1 ciclo da rede: N_AMOSTRAS amostras de tensão e de corrente,
// com o intervalo ajustado pela frequência medida (amostragem síncrona)
void coletarCiclo() {
  double intervaloUs = 1000000.0 / (freqRede * N_AMOSTRAS);
  unsigned long t0 = micros();
  for (int k = 0; k < N_AMOSTRAS; k++) {
    unsigned long alvo = t0 + (unsigned long)(k * intervaloUs);
    while ((long)(micros() - alvo) < 0) { }
    if ((long)(micros() - alvo) > (long)(intervaloUs / 2)) nAtrasos++;  // não deu tempo
    lerAmostra(amostrasV[k], amostrasI[k]);
  }
  duracaoBlocoUs = micros() - t0;
}

// FFT de um sinal: acumula o quadrado da magnitude de cada harmônico
// e devolve a componente complexa da fundamental (bin 1)
void analisarFFT(uint16_t *amostras, double *acumHarm2, double &re1, double &im1) {
  for (int k = 0; k < N_AMOSTRAS; k++) { vReal[k] = amostras[k]; vImag[k] = 0; }
  FFT.dcRemoval();                      // remove o offset
  FFT.compute(FFTDirection::Forward);   // sem janela: bloco com exatamente 1 ciclo
  re1 = vReal[1];
  im1 = vImag[1];
  FFT.complexToMagnitude();
  for (int h = 1; h <= N_HARMONICOS; h++) acumHarm2[h] += vReal[h] * vReal[h];
}

// Processa o ciclo coletado: RMS, potência, FFT e defasagem
void processarCiclo() {
  double mediaV = 0, mediaI = 0;
  for (int k = 0; k < N_AMOSTRAS; k++) { mediaV += amostrasV[k]; mediaI += amostrasI[k]; }
  mediaV /= N_AMOSTRAS;
  mediaI /= N_AMOSTRAS;
  mediaVAtual = mediaV;

  // Tensão sem offset, deslocada pela correção de fase (deslocamento circular,
  // válido porque o bloco tem exatamente 1 ciclo)
  double vCorr[N_AMOSTRAS];
  double desloc = CORRECAO_FASE_GRAUS / 360.0 * N_AMOSTRAS;
  for (int k = 0; k < N_AMOSTRAS; k++) {
    double pos = fmod(k - desloc, (double)N_AMOSTRAS);
    if (pos < 0) pos += N_AMOSTRAS;
    int k0 = (int)pos;
    int k1 = (k0 + 1) % N_AMOSTRAS;
    double fr = pos - k0;
    vCorr[k] = (amostrasV[k0] - mediaV) * (1 - fr) + (amostrasV[k1] - mediaV) * fr;
  }

  double sinalI = INVERTER_CORRENTE ? -1.0 : 1.0;
  for (int k = 0; k < N_AMOSTRAS; k++) {
    double cV = amostrasV[k] - mediaV;          // remove o offset
    double cI = amostrasI[k] - mediaI;
    double v = vCorr[k] * FATOR_TENSAO;          // tensão instantânea (V), com correção de fase
    double i = sinalI * cI * FATOR_CORRENTE / N_VOLTAS;  // corrente instantânea (A)
    somaBrutoV2 += cV * cV;
    somaBrutoI2 += cI * cI;
    somaV2 += v * v;
    somaI2 += i * i;
    somaP  += v * i;                             // potência instantânea
  }

  double vRe, vIm, iRe, iIm;
  analisarFFT(amostrasV, harmV2, vRe, vIm);
  analisarFFT(amostrasI, harmI2, iRe, iIm);
  if (INVERTER_CORRENTE) { iRe = -iRe; iIm = -iIm; }
  // V1 x conjugado(I1): o ângulo desse produto é a defasagem entre V e I
  cruzRe += vRe * iRe + vIm * iIm;
  cruzIm += vIm * iRe - vRe * iIm;

  nCiclos++;
}

// Mede a frequência pelos cruzamentos por zero (subida) da tensão
// numa janela de ~200 ms (12 ciclos), com interpolação linear
double medirFrequencia() {
  const unsigned long JANELA_US = 200000;
  const int HISTERESE = 60;            // contagens do ADC (evita falsos cruzamentos)
  double nivel = mediaVAtual;

  unsigned long inicio = micros();
  unsigned long tAnt = 0;
  int vAnt = lerADC(PINO_TENSAO);
  bool armado = false;
  double tPrimeira = 0, tUltima = 0;
  int nCruz = 0;

  while (true) {
    unsigned long t = micros() - inicio;
    if (t > JANELA_US) break;
    int v = lerADC(PINO_TENSAO);
    if (v < nivel - HISTERESE) armado = true;
    if (armado && vAnt < nivel && v >= nivel) {
      double frac = (nivel - vAnt) / (double)(v - vAnt);
      double tc = tAnt + frac * (double)(t - tAnt);
      if (nCruz == 0) tPrimeira = tc;
      tUltima = tc;
      nCruz++;
      armado = false;
    }
    vAnt = v;
    tAnt = t;
  }
  if (nCruz >= 3) return (nCruz - 1) * 1000000.0 / (tUltima - tPrimeira);
  return 0;   // sem sinal de tensão
}

// THD a partir dos valores RMS de cada harmônico (h[1] = fundamental)
double calcularTHD(double *h) {
  double soma = 0;
  for (int k = 2; k <= N_HARMONICOS; k++) soma += h[k] * h[k];
  return (h[1] > 0) ? sqrt(soma) / h[1] : 0;
}

// Calcula as grandezas do período e imprime o relatório
void imprimirRelatorio(double freqMedida) {
  double nTotal = (double)nCiclos * N_AMOSTRAS;

  double vrms = sqrt(somaV2 / nTotal);
  double irms = sqrt(somaI2 / nTotal);
  irms = sqrt(fmax(0.0, irms * irms - RUIDO_I_RMS * RUIDO_I_RMS));
  double p  = somaP / nTotal;
  double s  = vrms * irms;
  double q  = sqrt(fmax(0.0, s * s - p * p));
  double fp = (s > 0.001) ? p / s : 0;
  double defasagem = atan2(cruzIm, cruzRe) * 180.0 / PI - CORRECAO_FASE_GRAUS; // + indutivo
  while (defasagem > 180.0)  defasagem -= 360.0;
  while (defasagem < -180.0) defasagem += 360.0;
  double fpDesl = cos(defasagem * PI / 180.0);

  // Valor RMS de cada harmônico, convertido para volts e ampères
  // (magnitude da FFT -> RMS: multiplica por raiz(2) / N_AMOSTRAS)
  double hV[N_HARMONICOS + 1], hI[N_HARMONICOS + 1];
  for (int k = 1; k <= N_HARMONICOS; k++) {
    hV[k] = sqrt(harmV2[k] / nCiclos) * sqrt(2.0) / N_AMOSTRAS * FATOR_TENSAO;
    hI[k] = sqrt(harmI2[k] / nCiclos) * sqrt(2.0) / N_AMOSTRAS * FATOR_CORRENTE / N_VOLTAS;
  }
  double thdV = calcularTHD(hV);
  double thdI = calcularTHD(hI);
  double v1 = hV[1];   // tensão da fundamental (V)
  double i1 = hI[1];   // corrente da fundamental (A)

  // Há sinal quando a componente de 60 Hz supera o limiar
  // (o ruído se espalha por todas as frequências e quase não aparece no 60 Hz)
  bool temTensao   = v1 >= VRMS_MINIMO;
  bool temCorrente = i1 >= IRMS_MINIMO;

  Serial.println();
  Serial.println("==================== MEDICAO ====================");
  if (temTensao && freqMedida > 0) Serial.printf("Frequencia:          %.2f Hz\n", freqRede);
  else                             Serial.printf("Frequencia:          -- (sem tensao)\n");
  Serial.printf("Tensao RMS:          %.1f V   (fundamental %.1f V)\n", vrms, v1);
  Serial.printf("Corrente RMS:        %.3f A   (fundamental %.3f A)\n", irms, i1);

  if (temTensao && temCorrente) {
    Serial.printf("Potencia ativa P:    %.1f W\n", p);
    Serial.printf("Potencia reativa Q:  %.1f var\n", q);
    Serial.printf("Potencia aparente S: %.1f VA\n", s);
    const char *tipo = fabs(defasagem) < 1.0 ? "resistivo" : (defasagem > 0 ? "indutivo" : "capacitivo");
    Serial.printf("Fator de potencia:   %.3f (%s)\n", fp, tipo);
    Serial.printf("Defasagem V-I:       %.1f graus | cos(phi) = %.3f\n", defasagem, fpDesl);
    if (TRIFASICO_EQUILIBRADO) {
      Serial.printf("TOTAL TRIFASICO:     P = %.1f W | Q = %.1f var | S = %.1f VA\n", 3 * p, 3 * q, 3 * s);
    }
  } else {
    Serial.println("Potencias e fator de potencia: -- (sem tensao ou sem carga)");
  }

  if (temTensao) Serial.printf("THD tensao:          %.2f %%\n", thdV * 100);
  else           Serial.println("THD tensao:          --");
  if (temCorrente) Serial.printf("THD corrente:        %.2f %%\n", thdI * 100);
  else             Serial.println("THD corrente:        --");

  // Tabela da FFT
  Serial.println("---------- FFT: harmonicos (% da fundamental) ----------");
  Serial.println(" H   Freq     Tensao (V)   %V      Corrente (A)   %I");
  for (int k = 1; k <= N_HARMONICOS; k++) {
    Serial.printf("%2d  %4.0f Hz ", k, k * freqRede);
    if (temTensao)   Serial.printf("  %8.2f   %6.2f  ", hV[k], hV[k] / hV[1] * 100);
    else             Serial.printf("        --       --  ");
    if (temCorrente) Serial.printf("  %8.3f     %6.2f\n", hI[k], hI[k] / hI[1] * 100);
    else             Serial.printf("        --         --\n");
  }

  Serial.printf("[calibracao] RMS bruto V: %.1f | RMS bruto I: %.1f | ciclos: %d\n",
                sqrt(somaBrutoV2 / nTotal), sqrt(somaBrutoI2 / nTotal), nCiclos);
  Serial.printf("[amostragem] bloco de %.2f ms (ideal %.2f ms) | amostras atrasadas: %d\n",
                duracaoBlocoUs / 1000.0, 1000.0 / freqRede, nAtrasos);
  if (nAtrasos > 0) {
    Serial.println("AVISO: o ADC nao acompanhou a amostragem. Reduza N_AMOSTRAS para 32.");
    Serial.println("       Os harmonicos e a THD deste relatorio nao sao confiaveis.");
  }
}

// ================================================================
void setup() {
  Serial.begin(115200);
  analogReadResolution(12);         // ADC de 12 bits (0 a 4095)
  analogSetAttenuation(ADC_11db);   // faixa de 0 a ~3,3 V
  delay(1000);
  Serial.println("Prototipo de monitoramento de qualidade de energia iniciado.");

  // Mede quanto tempo o ADC leva para ler uma amostra de tensão + corrente
  unsigned long t = micros();
  uint16_t va, ia;
  for (int n = 0; n < 100; n++) lerAmostra(va, ia);
  double tAmostra = (micros() - t) / 100.0;
  double tDisponivel = 1000000.0 / (60.0 * N_AMOSTRAS);
  Serial.printf("Tempo por amostra (V + I): %.1f us | intervalo disponivel: %.1f us\n",
                tAmostra, tDisponivel);
  if (tAmostra > 0.9 * tDisponivel)
    Serial.println("AVISO: ADC lento para esta configuracao. Reduza N_AMOSTRAS para 32.");
}

void loop() {
  // 1) Mede N_CICLOS_MEDIA ciclos de tensão e corrente
  zerarAcumuladores();
  for (int c = 0; c < N_CICLOS_MEDIA; c++) {
    coletarCiclo();
    processarCiclo();
  }

  // 2) Mede a frequência e ajusta o sincronismo da amostragem
  double f = medirFrequencia();
  if (f > 55.0 && f < 65.0) freqRede = 0.7 * freqRede + 0.3 * f;

  // 3) Mostra os resultados
  imprimirRelatorio(f);
}
