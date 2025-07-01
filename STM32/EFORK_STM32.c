/* USER CODE BEGIN Header */
/**
 Archivo: main.c - EFORK (orden fraccionario) en STM32
 Salida a DAC (canales 1 y 2) y envío UART asíncrono (interrupción).
**/

/* USER CODE BEGIN Header */
/**
 *   ******************************************************************************
  * @file           : main.c
  * @brief          : EFORK para Lorenz (orden fraccionario) en STM32F7,
  *                   salidas DAC1/DAC2 y UART asíncrona (IT).
  ******************************************************************************
  */
/* USER CODE END Header */
/* Includes ------------------------------------------------------------------*/
#include "main.h"          
#include <math.h>
#include <string.h>
#include <stdio.h>
/* Private includes ----------------------------------------------------------*/
/* USER CODE BEGIN Includes */

/* USER CODE END Includes */

/* Private typedef -----------------------------------------------------------*/
/* USER CODE BEGIN PTD */

/* USER CODE END PTD */

/* Private define ------------------------------------------------------------*/
/* USER CODE BEGIN PD */
// Ajusta según tu orden fraccionario y longitud de memoria deseada
#define ALPHA   0.995
#define DT      0.005
#define Lm      10.0

// Parámetros Lorenz
#define SIGMA     10.0
#define RHO       28.0
#define BETA      (8.0/3.0)

// Máximo tamaño de ventana de memoria
#define MAX_M     10000  // Ajusta según RAM disponible

/* Parámetros de los tres sistemas (Lorenz, Rossler, Chen) */
/* #define SIGMA   10.0
#define RHO     28.0
#define BETA    (8.0/3.0)
#define U       7.5
#define V       1.0
#define W       5.0
#define A       0.2
#define BB      0.2      //“B” ya reservada por HAL 
#define CC      5.7 
*/
// Dimensión del sistema Lorenz
#define D         3
/* USER CODE END PD */

/* Private macro -------------------------------------------------------------*/
/* USER CODE BEGIN PM */

/* USER CODE END PM */

/* Private variables ---------------------------------------------------------*/
#if defined ( __ICCARM__ ) /*!< IAR Compiler */
#pragma location=0x2004c000
ETH_DMADescTypeDef  DMARxDscrTab[ETH_RX_DESC_CNT]; /* Ethernet Rx DMA Descriptors */
#pragma location=0x2004c0a0
ETH_DMADescTypeDef  DMATxDscrTab[ETH_TX_DESC_CNT]; /* Ethernet Tx DMA Descriptors */

#elif defined ( __CC_ARM )  /* MDK ARM Compiler */

__attribute__((at(0x2004c000))) ETH_DMADescTypeDef  DMARxDscrTab[ETH_RX_DESC_CNT]; /* Ethernet Rx DMA Descriptors */
__attribute__((at(0x2004c0a0))) ETH_DMADescTypeDef  DMATxDscrTab[ETH_TX_DESC_CNT]; /* Ethernet Tx DMA Descriptors */

#elif defined ( __GNUC__ ) /* GNU Compiler */

ETH_DMADescTypeDef DMARxDscrTab[ETH_RX_DESC_CNT] __attribute__((section(".RxDecripSection"))); /* Ethernet Rx DMA Descriptors */
ETH_DMADescTypeDef DMATxDscrTab[ETH_TX_DESC_CNT] __attribute__((section(".TxDecripSection")));   /* Ethernet Tx DMA Descriptors */
#endif

ETH_TxPacketConfig TxConfig;

DAC_HandleTypeDef hdac;
DMA_HandleTypeDef hdma_dac1;
DMA_HandleTypeDef hdma_dac2;

ETH_HandleTypeDef heth;

I2C_HandleTypeDef hi2c1;

TIM_HandleTypeDef htim2;

UART_HandleTypeDef huart3;
DMA_HandleTypeDef hdma_usart3_rx;
DMA_HandleTypeDef hdma_usart3_tx;

PCD_HandleTypeDef hpcd_USB_OTG_FS;
/* USER CODE BEGIN PV */

// Longitud de memoria (en segundos), y número de muestras m = Lm / DT
static int    m = 0;      // Se calculará en el main

// Coeficientes binomiales c[j], 0 <= j <= m
//static long double c_coef[MAX_M + 1];

// Historial de las variables x, y, z (memoria corta)
static long double x_hist[MAX_M + 1];
static long double y_hist[MAX_M + 1];
static long double z_hist[MAX_M + 1];

static long double t_hist[MAX_M+1];                // NUEVO (para EFORK)

// Índice circular para indicar la posición más reciente en el historial
static int idx_current = 0;

// Variables para enviar por UART sin bloqueo
static char tx_buffer[100];
static volatile uint8_t uart_busy = 0; // 0=libre, 1=ocupada

// Contador de pasos de integración
static int step_count = 0;

/* USER CODE END PV */

/* Private function prototypes -----------------------------------------------*/
void SystemClock_Config(void);
static void MX_GPIO_Init(void);
static void MX_DMA_Init(void);
static void MX_ETH_Init(void);
static void MX_I2C1_Init(void);
static void MX_USART3_UART_Init(void);
static void MX_USB_OTG_FS_PCD_Init(void);
static void MX_DAC_Init(void);
static void MX_TIM2_Init(void);
/* USER CODE BEGIN PFP */

// 1) Función para calcular coeficientes binomiales para GL
//static void gl_compute_binomial_coefs(long double alpha, int m);

// 2) Paso de integración del sistema de Lorenz usando EFORK (memoria corta)
//static void lorenz_efork_step(void);

// 3) Conversión de float/double a DAC 12 bits (rango [-20..20])
static uint32_t float_to_dac12(long double val, long double minVal, long double maxVal);

// 4) Envío no bloqueante por UART
static void enviar_mensaje(const char *msg);

// 5) Callback de final de envío UART
void HAL_UART_TxCpltCallback(UART_HandleTypeDef *huart);

/* USER CODE END PFP */
/* Private user code ---------------------------------------------------------*/
// ==============================================================================
// No se usa para EFORK, modificar para lo que se necesite calcular para este nuevo método.
/* USER CODE BEGIN 0 */
/**
  * @brief  Calcula los coeficientes binomiales c[j] para la aproximación
  *         de Grünwald–Letnikov con una relación recurrente.
  * @note   Solo se hace una vez al inicio, para j=0..m.
  */

//static void gl_compute_binomial_coefs(long double alpha, int m)
/* {
    c_coef[0] = 1.0;
    for(int j = 1; j <= m; j++)
    {
        // c[j] = c[j-1] - (c[j-1]*(1+alpha)/j)
    	long double tmp = (1.0 + alpha) / (long double)j;
    	//c_coef[j] = c_coef[j - 1] - (c_coef[j - 1] * tmp);
    	c_coef[j] = c_coef[j - 1]*(1 - tmp) ;

    }
} */
//==================================================================================
/***** Memoria fraccional (ventana corta m) *****/
static long double memory_fractional(int k,long double t,
                                     const long double *v,const long double *tt)
{
    if(k==0) return 0.0L;
    long double gamma_term = tgammal(2.0L - ALPHA);
    int start = (k-m>0)? (k-m):0;
    long double sum = 0.0L;
    for(int j=start;j<k;j++){
        long double vdiff = v[j+1]-v[j];
        long double t0 = tt[j];
        long double t1 = tt[j+1];
        long double f1 = powl(t - t0, 1.0L-ALPHA);
        long double f2 = powl(t - t1, 1.0L-ALPHA);
        sum += vdiff*(f1-f2);
    }
    return sum/(DT*gamma_term);
}

/***** Paso EFORK (RK3‑α) *****/
static void efork_step(void)
{
    /* --------- calcular constantes EFORK una sola vez --------- */
    static int init = 0;
    static long double ha,g1,g2,g3,a21,a31,a32,w1,w2,w3;
    if(!init){
        g1 = tgammal(1.0L+ALPHA);
        g2 = tgammal(1.0L+2.0L*ALPHA);
        g3 = tgammal(1.0L+3.0L*ALPHA);
        ha = powl(DT,ALPHA);
        a21 = 1.0L/(2.0L*g1*g1);
        a31 = (g1*g1*tgammal(2*ALPHA+1)+2.0L*tgammal(2*ALPHA+1)*tgammal(2*ALPHA+1)-tgammal(3*ALPHA+1)) /
               (4.0L*g1*g1*(2.0L*tgammal(2*ALPHA+1)*tgammal(2*ALPHA+1)-tgammal(3*ALPHA+1)));
        a32 = -tgammal(2*ALPHA+1)/(4.0L*(2.0L*tgammal(2*ALPHA+1)*tgammal(2*ALPHA+1)-tgammal(3*ALPHA+1)));
        w1 = (8.0L*powl(g1,3)*powl(tgammal(1+2*ALPHA),2) - 6.0L*powl(g1,3)*tgammal(1+3*ALPHA) + tgammal(1+2*ALPHA)*tgammal(1+3*ALPHA)) /
              (g1*tgammal(1+2*ALPHA)*tgammal(1+3*ALPHA));
        w2 = 2.0L*g1*g1*(4.0L*powl(tgammal(1+2*ALPHA),2) - tgammal(1+3*ALPHA)) /
              (tgammal(1+2*ALPHA)*tgammal(1+3*ALPHA));
        w3 = -8.0L*g1*g1*(2.0L*powl(tgammal(1+2*ALPHA),2) - tgammal(1+3*ALPHA)) /
               (tgammal(1+2*ALPHA)*tgammal(1+3*ALPHA));
        init = 1;
    }

    /* ------ estado actual ------ */
    long double t_n = t_hist[idx_current];
    long double x_n = x_hist[idx_current];
    long double y_n = y_hist[idx_current];
    long double z_n = z_hist[idx_current];

    int k = step_count;
    long double mem_x = memory_fractional(k,t_n,x_hist,t_hist);
    long double mem_y = memory_fractional(k,t_n,y_hist,t_hist);
    long double mem_z = memory_fractional(k,t_n,z_hist,t_hist);

    long double dx,dy,dz;
    /* K1 */
    system_rhs(x_n,y_n,z_n,&dx,&dy,&dz);
    long double K1x = ha*(dx - mem_x);
    long double K1y = ha*(dy - mem_y);
    long double K1z = ha*(dz - mem_z);
    /* K2 */
    system_rhs(x_n + a21*K1x, y_n + a21*K1y, z_n + a21*K1z, &dx,&dy,&dz);
    long double K2x = ha*dx;
    long double K2y = ha*dy;
    long double K2z = ha*dz;
    /* K3 */
    system_rhs(x_n + a31*K2x + a32*K1x,
               y_n + a31*K2y + a32*K1y,
               z_n + a31*K2z + a32*K1z, &dx,&dy,&dz);
    long double K3x = ha*dx;
    long double K3y = ha*dy;
    long double K3z = ha*dz;

    /* nuevo estado */
    long double x_np1 = x_n + w1*K1x + w2*K2x + w3*K3x;
    long double y_np1 = y_n + w1*K1y + w2*K2y + w3*K3y;
    long double z_np1 = z_n + w1*K1z + w2*K2z + w3*K3z;

    /* guardar en buffer */
    int next = (idx_current+1)%(MAX_M+1);
    t_hist[next] = t_n + DT;
    x_hist[next] = x_np1;
    y_hist[next] = y_np1;
    z_hist[next] = z_np1;
    idx_current  = next;
    step_count++;

    /* salidas */
    HAL_DAC_SetValue(&hdac,DAC_CHANNEL_1,DAC_ALIGN_12B_R,float_to_dac12(x_np1,-25.0L,25.0L));
    HAL_DAC_SetValue(&hdac,DAC_CHANNEL_2,DAC_ALIGN_12B_R,float_to_dac12(y_np1,-30.0L,30.0L));
    char buf[90];
    snprintf(buf,sizeof(buf),"EFORK -> x=% .3Lf, y=% .3Lf, z=% .3Lf\r\n",x_np1,y_np1,z_np1);
    enviar_msj(buf);
}

/***** UART helper (igual que en GL salvo nombre función) *****/
static void enviar_msj(const char *msg) { /* = copia igual que en GL, solo cambia nombre */ }

void HAL_UART_TxCpltCallback(UART_HandleTypeDef *huart){ /* = igual que en GL */ }

/***** float_to_dac12 = igual que en GL *****/

/* ---------- 6. main() : SOLO las líneas distintas ---------- */
int main(void)
{
    HAL_Init();
    SystemClock_Config();       // = mismo que en GL
    MX_GPIO_Init();             // = mismo que en GL
    MX_DMA_Init();              // = mismo que en GL
    MX_DAC_Init();              // = mismo que en GL
    MX_USART3_UART_Init();      // = mismo que en GL

    /* ---- calculamos m y llenamos históricos ---- */
    m = (int)(LM/DT); if(m>MAX_M) m = MAX_M;
    for(int i=0;i<=m;i++){ t_hist[i] = -DT*i; x_hist[i]=y_hist[i]=z_hist[i]=0.0L; }
    x_hist[0]=0.1L; y_hist[0]=0.1L; z_hist[0]=0.1L; t_hist[0]=0.0L;
    idx_current=0; step_count=0;

    HAL_DAC_Start(&hdac,DAC_CHANNEL_1);
    HAL_DAC_Start(&hdac,DAC_CHANNEL_2);
    enviar_msj("\r\n*** EFORK Frac + UART IT ***\r\n");

    while(1){
        efork_step();
        HAL_Delay(10);   /* ~100 Hz */
    }
}

/* ---------- El resto de funciones (SystemClock_Config, MX_xxx_Init, etc.)
 * ------------ son EXACTAMENTE las mismas que en tu archivo GL -------------- */

/**
  * @brief  The application entry point.
  * @retval int
  */
int main(void)
{

  /* USER CODE BEGIN 1 */

  /* USER CODE END 1 */

  /* MCU Configuration--------------------------------------------------------*/

  /* Reset of all peripherals, Initializes the Flash interface and the Systick. */
  HAL_Init();

  /* USER CODE BEGIN Init */

  /* USER CODE END Init */

  /* Configure the system clock */
  SystemClock_Config();

  /* USER CODE BEGIN SysInit */

  /* USER CODE END SysInit */

  /* Initialize all configured peripherals */
  MX_GPIO_Init();
  MX_DMA_Init();
  MX_ETH_Init();
  MX_I2C1_Init();
  MX_USART3_UART_Init();
  MX_USB_OTG_FS_PCD_Init();
  MX_DAC_Init();
  MX_TIM2_Init();

  /* USER CODE BEGIN 2 */
  HAL_TIM_Base_Start(&htim2);

    /* USER CODE BEGIN 2 */

    // 1) Calcula m = Lm / DT, y si excede MAX_M lo limitamos
    m = (int)(Lm / DT);
    if(m > MAX_M) m = MAX_M;

    // 2) Calculamos los coeficientes binomiales
    gl_compute_binomial_coefs(ALPHA, m);


    // 3) Inicializamos el historial COMPLETO en cero
    for(int i = 0; i <= m; i++)
    {
        x_hist[i] = 0.0;
        y_hist[i] = 0.0;
        z_hist[i] = 0.0;
    }


    // Luego pones tu condición inicial sólo en [0]
    long double x0 = 0.1, y0 = 0.1, z0 = 0.1;
    x_hist[0] = x0;
    y_hist[0] = y0;
    z_hist[0] = z0;

    // Arrancamos con idx_current = 0
    idx_current = 0;

    // Contador de pasos, empieza en 0
    step_count = 0;

    // Arrancamos DAC (canales 1 y 2)
    HAL_DAC_Start(&hdac, DAC_CHANNEL_1);
    HAL_DAC_Start(&hdac, DAC_CHANNEL_2);

    // Mensaje inicial por UART
    enviar_mensaje("\r\n*** Lorenz Frac (GL) + UART IT ***\r\n");

    // Opción: iniciar un timer base si se desea (ej. htim2)
    HAL_TIM_Base_Start(&htim2);

  /* USER CODE END 2 */

  /* Infinite loop */
  /* USER CODE BEGIN WHILE */
    while (1)
    {
    	// 1) Un paso de integración (GL)
    		      lorenz_gl_step();

    		      // 2) Obtener valores actuales
    		      long double x = x_hist[idx_current];
    		      long double y = y_hist[idx_current];
    		      long double z = z_hist[idx_current];

    		      // 3) Actualiza DAC
    		      uint32_t dac_x = float_to_dac12(x, -20.0L, 20.0L);
    		      uint32_t dac_y = float_to_dac12(y, -25.0L, 25.0L);
    		      //uint32_t dac_z = float_to_dac12(z, 15.0L, 50.0L);


    		      HAL_DAC_SetValue(&hdac, DAC_CHANNEL_1, DAC_ALIGN_12B_R, dac_x);
    		      HAL_DAC_SetValue(&hdac, DAC_CHANNEL_2, DAC_ALIGN_12B_R, dac_y);


    		      // 4) Envía por UART (sin bloquear)
    		      char buf[60];
    		      snprintf(buf, sizeof(buf),
    		               "FracGL -> x=%.3Lf, y=%.3Lf, z=%.3Lf\r\n", x, y, z);
    		      enviar_mensaje(buf);

    		      // 5) Pequeño retardo => ~100 Hz
    		      HAL_Delay(10);

    }
  /* USER CODE END 3 */
}

/**
  * @brief System Clock Configuration
  * @retval None
  */
void SystemClock_Config(void)
{
  RCC_OscInitTypeDef RCC_OscInitStruct = {0};
  RCC_ClkInitTypeDef RCC_ClkInitStruct = {0};

  /** Configure LSE Drive Capability
  */
  HAL_PWR_EnableBkUpAccess();

  /** Configure the main internal regulator output voltage
  */
  __HAL_RCC_PWR_CLK_ENABLE();
  __HAL_PWR_VOLTAGESCALING_CONFIG(PWR_REGULATOR_VOLTAGE_SCALE3);

  /** Initializes the RCC Oscillators according to the specified parameters
  * in the RCC_OscInitTypeDef structure.
  */
  RCC_OscInitStruct.OscillatorType = RCC_OSCILLATORTYPE_HSE;
  RCC_OscInitStruct.HSEState = RCC_HSE_BYPASS;
  RCC_OscInitStruct.PLL.PLLState = RCC_PLL_ON;
  RCC_OscInitStruct.PLL.PLLSource = RCC_PLLSOURCE_HSE;
  RCC_OscInitStruct.PLL.PLLM = 4;
  RCC_OscInitStruct.PLL.PLLN = 72;
  RCC_OscInitStruct.PLL.PLLP = RCC_PLLP_DIV2;
  RCC_OscInitStruct.PLL.PLLQ = 3;
  if (HAL_RCC_OscConfig(&RCC_OscInitStruct) != HAL_OK)
  {
    Error_Handler();
  }

  /** Initializes the CPU, AHB and APB buses clocks
  */
  RCC_ClkInitStruct.ClockType = RCC_CLOCKTYPE_HCLK|RCC_CLOCKTYPE_SYSCLK
                              |RCC_CLOCKTYPE_PCLK1|RCC_CLOCKTYPE_PCLK2;
  RCC_ClkInitStruct.SYSCLKSource = RCC_SYSCLKSOURCE_PLLCLK;
  RCC_ClkInitStruct.AHBCLKDivider = RCC_SYSCLK_DIV1;
  RCC_ClkInitStruct.APB1CLKDivider = RCC_HCLK_DIV2;
  RCC_ClkInitStruct.APB2CLKDivider = RCC_HCLK_DIV1;

  if (HAL_RCC_ClockConfig(&RCC_ClkInitStruct, FLASH_LATENCY_2) != HAL_OK)
  {
    Error_Handler();
  }
}

/**
  * @brief DAC Initialization Function
  * @param None
  * @retval None
  */
static void MX_DAC_Init(void)
{

  /* USER CODE BEGIN DAC_Init 0 */

  /* USER CODE END DAC_Init 0 */

  DAC_ChannelConfTypeDef sConfig = {0};

  /* USER CODE BEGIN DAC_Init 1 */

  /* USER CODE END DAC_Init 1 */

  /** DAC Initialization
  */
  hdac.Instance = DAC;
  if (HAL_DAC_Init(&hdac) != HAL_OK)
  {
    Error_Handler();
  }

  /** DAC channel OUT1 config
  */
  sConfig.DAC_Trigger = DAC_TRIGGER_T2_TRGO;
  sConfig.DAC_OutputBuffer = DAC_OUTPUTBUFFER_ENABLE;
  if (HAL_DAC_ConfigChannel(&hdac, &sConfig, DAC_CHANNEL_1) != HAL_OK)
  {
    Error_Handler();
  }

  /** DAC channel OUT2 config
  */
  if (HAL_DAC_ConfigChannel(&hdac, &sConfig, DAC_CHANNEL_2) != HAL_OK)
  {
    Error_Handler();
  }
  /* USER CODE BEGIN DAC_Init 2 */

  /* USER CODE END DAC_Init 2 */

}

/**
  * @brief ETH Initialization Function
  * @param None
  * @retval None
  */
static void MX_ETH_Init(void)
{

  /* USER CODE BEGIN ETH_Init 0 */

  /* USER CODE END ETH_Init 0 */

   static uint8_t MACAddr[6];

  /* USER CODE BEGIN ETH_Init 1 */

  /* USER CODE END ETH_Init 1 */
  heth.Instance = ETH;
  MACAddr[0] = 0x00;
  MACAddr[1] = 0x80;
  MACAddr[2] = 0xE1;
  MACAddr[3] = 0x00;
  MACAddr[4] = 0x00;
  MACAddr[5] = 0x00;
  heth.Init.MACAddr = &MACAddr[0];
  heth.Init.MediaInterface = HAL_ETH_RMII_MODE;
  heth.Init.TxDesc = DMATxDscrTab;
  heth.Init.RxDesc = DMARxDscrTab;
  heth.Init.RxBuffLen = 1524;

  /* USER CODE BEGIN MACADDRESS */

  /* USER CODE END MACADDRESS */

  if (HAL_ETH_Init(&heth) != HAL_OK)
  {
    Error_Handler();
  }

  memset(&TxConfig, 0 , sizeof(ETH_TxPacketConfig));
  TxConfig.Attributes = ETH_TX_PACKETS_FEATURES_CSUM | ETH_TX_PACKETS_FEATURES_CRCPAD;
  TxConfig.ChecksumCtrl = ETH_CHECKSUM_IPHDR_PAYLOAD_INSERT_PHDR_CALC;
  TxConfig.CRCPadCtrl = ETH_CRC_PAD_INSERT;
  /* USER CODE BEGIN ETH_Init 2 */

  /* USER CODE END ETH_Init 2 */

}

/**
  * @brief I2C1 Initialization Function
  * @param None
  * @retval None
  */
static void MX_I2C1_Init(void)
{

  /* USER CODE BEGIN I2C1_Init 0 */

  /* USER CODE END I2C1_Init 0 */

  /* USER CODE BEGIN I2C1_Init 1 */

  /* USER CODE END I2C1_Init 1 */
  hi2c1.Instance = I2C1;
  hi2c1.Init.Timing = 0x00808CD2;
  hi2c1.Init.OwnAddress1 = 0;
  hi2c1.Init.AddressingMode = I2C_ADDRESSINGMODE_7BIT;
  hi2c1.Init.DualAddressMode = I2C_DUALADDRESS_DISABLE;
  hi2c1.Init.OwnAddress2 = 0;
  hi2c1.Init.OwnAddress2Masks = I2C_OA2_NOMASK;
  hi2c1.Init.GeneralCallMode = I2C_GENERALCALL_DISABLE;
  hi2c1.Init.NoStretchMode = I2C_NOSTRETCH_DISABLE;
  if (HAL_I2C_Init(&hi2c1) != HAL_OK)
  {
    Error_Handler();
  }

  /** Configure Analogue filter
  */
  if (HAL_I2CEx_ConfigAnalogFilter(&hi2c1, I2C_ANALOGFILTER_ENABLE) != HAL_OK)
  {
    Error_Handler();
  }

  /** Configure Digital filter
  */
  if (HAL_I2CEx_ConfigDigitalFilter(&hi2c1, 0) != HAL_OK)
  {
    Error_Handler();
  }
  /* USER CODE BEGIN I2C1_Init 2 */

  /* USER CODE END I2C1_Init 2 */

}

/**
  * @brief TIM2 Initialization Function
  * @param None
  * @retval None
  */
static void MX_TIM2_Init(void)
{

  /* USER CODE BEGIN TIM2_Init 0 */

  /* USER CODE END TIM2_Init 0 */

  TIM_ClockConfigTypeDef sClockSourceConfig = {0};
  TIM_MasterConfigTypeDef sMasterConfig = {0};

  /* USER CODE BEGIN TIM2_Init 1 */

  /* USER CODE END TIM2_Init 1 */
  htim2.Instance = TIM2;
  htim2.Init.Prescaler = 72-1;
  htim2.Init.CounterMode = TIM_COUNTERMODE_UP;
  htim2.Init.Period = 100;
  htim2.Init.ClockDivision = TIM_CLOCKDIVISION_DIV1;
  htim2.Init.AutoReloadPreload = TIM_AUTORELOAD_PRELOAD_DISABLE;
  if (HAL_TIM_Base_Init(&htim2) != HAL_OK)
  {
    Error_Handler();
  }
  sClockSourceConfig.ClockSource = TIM_CLOCKSOURCE_INTERNAL;
  if (HAL_TIM_ConfigClockSource(&htim2, &sClockSourceConfig) != HAL_OK)
  {
    Error_Handler();
  }
  sMasterConfig.MasterOutputTrigger = TIM_TRGO_UPDATE;
  sMasterConfig.MasterSlaveMode = TIM_MASTERSLAVEMODE_DISABLE;
  if (HAL_TIMEx_MasterConfigSynchronization(&htim2, &sMasterConfig) != HAL_OK)
  {
    Error_Handler();
  }
  /* USER CODE BEGIN TIM2_Init 2 */

  /* USER CODE END TIM2_Init 2 */

}

/**
  * @brief USART3 Initialization Function
  * @param None
  * @retval None
  */
static void MX_USART3_UART_Init(void)
{

  /* USER CODE BEGIN USART3_Init 0 */

  /* USER CODE END USART3_Init 0 */

  /* USER CODE BEGIN USART3_Init 1 */

  /* USER CODE END USART3_Init 1 */
  huart3.Instance = USART3;
  huart3.Init.BaudRate = 115200;
  huart3.Init.WordLength = UART_WORDLENGTH_8B;
  huart3.Init.StopBits = UART_STOPBITS_1;
  huart3.Init.Parity = UART_PARITY_NONE;
  huart3.Init.Mode = UART_MODE_TX_RX;
  huart3.Init.HwFlowCtl = UART_HWCONTROL_NONE;
  huart3.Init.OverSampling = UART_OVERSAMPLING_16;
  huart3.Init.OneBitSampling = UART_ONE_BIT_SAMPLE_DISABLE;
  huart3.AdvancedInit.AdvFeatureInit = UART_ADVFEATURE_NO_INIT;
  if (HAL_UART_Init(&huart3) != HAL_OK)
  {
    Error_Handler();
  }
  /* USER CODE BEGIN USART3_Init 2 */

  /* USER CODE END USART3_Init 2 */

}

/**
  * @brief USB_OTG_FS Initialization Function
  * @param None
  * @retval None
  */
static void MX_USB_OTG_FS_PCD_Init(void)
{

  /* USER CODE BEGIN USB_OTG_FS_Init 0 */

  /* USER CODE END USB_OTG_FS_Init 0 */

  /* USER CODE BEGIN USB_OTG_FS_Init 1 */

  /* USER CODE END USB_OTG_FS_Init 1 */
  hpcd_USB_OTG_FS.Instance = USB_OTG_FS;
  hpcd_USB_OTG_FS.Init.dev_endpoints = 6;
  hpcd_USB_OTG_FS.Init.speed = PCD_SPEED_FULL;
  hpcd_USB_OTG_FS.Init.dma_enable = DISABLE;
  hpcd_USB_OTG_FS.Init.phy_itface = PCD_PHY_EMBEDDED;
  hpcd_USB_OTG_FS.Init.Sof_enable = ENABLE;
  hpcd_USB_OTG_FS.Init.low_power_enable = DISABLE;
  hpcd_USB_OTG_FS.Init.lpm_enable = DISABLE;
  hpcd_USB_OTG_FS.Init.vbus_sensing_enable = ENABLE;
  hpcd_USB_OTG_FS.Init.use_dedicated_ep1 = DISABLE;
  if (HAL_PCD_Init(&hpcd_USB_OTG_FS) != HAL_OK)
  {
    Error_Handler();
  }
  /* USER CODE BEGIN USB_OTG_FS_Init 2 */

  /* USER CODE END USB_OTG_FS_Init 2 */

}

/**
  * Enable DMA controller clock
  */
static void MX_DMA_Init(void)
{

  /* DMA controller clock enable */
  __HAL_RCC_DMA1_CLK_ENABLE();

  /* DMA interrupt init */
  /* DMA1_Stream1_IRQn interrupt configuration */
  HAL_NVIC_SetPriority(DMA1_Stream1_IRQn, 0, 0);
  HAL_NVIC_EnableIRQ(DMA1_Stream1_IRQn);
  /* DMA1_Stream3_IRQn interrupt configuration */
  HAL_NVIC_SetPriority(DMA1_Stream3_IRQn, 0, 0);
  HAL_NVIC_EnableIRQ(DMA1_Stream3_IRQn);
  /* DMA1_Stream5_IRQn interrupt configuration */
  HAL_NVIC_SetPriority(DMA1_Stream5_IRQn, 0, 0);
  HAL_NVIC_EnableIRQ(DMA1_Stream5_IRQn);
  /* DMA1_Stream6_IRQn interrupt configuration */
  HAL_NVIC_SetPriority(DMA1_Stream6_IRQn, 0, 0);
  HAL_NVIC_EnableIRQ(DMA1_Stream6_IRQn);

}

/**
  * @brief GPIO Initialization Function
  * @param None
  * @retval None
  */
static void MX_GPIO_Init(void)
{
  GPIO_InitTypeDef GPIO_InitStruct = {0};
  /* USER CODE BEGIN MX_GPIO_Init_1 */

  /* USER CODE END MX_GPIO_Init_1 */

  /* GPIO Ports Clock Enable */
  __HAL_RCC_GPIOC_CLK_ENABLE();
  __HAL_RCC_GPIOH_CLK_ENABLE();
  __HAL_RCC_GPIOA_CLK_ENABLE();
  __HAL_RCC_GPIOB_CLK_ENABLE();
  __HAL_RCC_GPIOD_CLK_ENABLE();
  __HAL_RCC_GPIOG_CLK_ENABLE();

  /*Configure GPIO pin Output Level */
  HAL_GPIO_WritePin(GPIOB, LD1_Pin|LD3_Pin|LD2_Pin, GPIO_PIN_RESET);

  /*Configure GPIO pin Output Level */
  HAL_GPIO_WritePin(USB_PowerSwitchOn_GPIO_Port, USB_PowerSwitchOn_Pin, GPIO_PIN_RESET);

  /*Configure GPIO pin : USER_Btn_Pin */
  GPIO_InitStruct.Pin = USER_Btn_Pin;
  GPIO_InitStruct.Mode = GPIO_MODE_IT_RISING;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  HAL_GPIO_Init(USER_Btn_GPIO_Port, &GPIO_InitStruct);

  /*Configure GPIO pins : LD1_Pin LD3_Pin LD2_Pin */
  GPIO_InitStruct.Pin = LD1_Pin|LD3_Pin|LD2_Pin;
  GPIO_InitStruct.Mode = GPIO_MODE_OUTPUT_PP;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  GPIO_InitStruct.Speed = GPIO_SPEED_FREQ_LOW;
  HAL_GPIO_Init(GPIOB, &GPIO_InitStruct);

  /*Configure GPIO pin : USB_PowerSwitchOn_Pin */
  GPIO_InitStruct.Pin = USB_PowerSwitchOn_Pin;
  GPIO_InitStruct.Mode = GPIO_MODE_OUTPUT_PP;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  GPIO_InitStruct.Speed = GPIO_SPEED_FREQ_LOW;
  HAL_GPIO_Init(USB_PowerSwitchOn_GPIO_Port, &GPIO_InitStruct);

  /*Configure GPIO pin : USB_OverCurrent_Pin */
  GPIO_InitStruct.Pin = USB_OverCurrent_Pin;
  GPIO_InitStruct.Mode = GPIO_MODE_INPUT;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  HAL_GPIO_Init(USB_OverCurrent_GPIO_Port, &GPIO_InitStruct);

  /* USER CODE BEGIN MX_GPIO_Init_2 */

  /* USER CODE END MX_GPIO_Init_2 */
}

/* USER CODE BEGIN 4 */

/* USER CODE END 4 */

/**
  * @brief  This function is executed in case of error occurrence.
  * @retval None
  */
void Error_Handler(void)
{
  /* USER CODE BEGIN Error_Handler_Debug */
  __disable_irq();
  while (1)
  {
  }
  /* USER CODE END Error_Handler_Debug */
}

#ifdef  USE_FULL_ASSERT
/**
  * @brief  Reports the name of the source file and the source line number
  *         where the assert_param error has occurred.
  * @param  file: pointer to the source file name
  * @param  line: assert_param error line source number
  * @retval None
  */
void assert_failed(uint8_t *file, uint32_t line)
{
  /* USER CODE BEGIN 6 */
  /* User can add his own implementation to report the file name and line number,
     ex: printf("Wrong parameters value: file %s on line %d\r\n", file, line) */
  /* USER CODE END 6 */
}
#endif /* USE_FULL_ASSERT */
