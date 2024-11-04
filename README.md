# BASE: a web service for providing compound-protein binding affinity prediction datasets with reduced similarity bias

## Overview
Deep learning-based drug-target affinity (DTA) prediction models have shown high performance but suffer from dataset bias. Our study investigates this bias using comprehensive databases and demonstrates that compound-protein binding affinity can often be predicted using compound features alone, due to high similarity among target proteins. We developed bias-reduced datasets by decreasing protein similarity between training and test sets, which improved model performance and balanced feature importance.

![figure](https://github.com/user-attachments/assets/96bf2278-4e52-4f4d-bf3d-f2b97484a0b2)

We introduce the Binding Affinity Similarity Explorer (BASE) web service, which offers bias-reduced datasets and prediction results to aid in the development of generalized and robust DTA models. BASE is freely available at https://synbi2024.kaist.ac.kr/base.

![Figure9](https://github.com/user-attachments/assets/18336449-4f46-4ffd-9164-b80008e92d5a)

## Installation Instructions
To run the project locally, clone the repository:
   ```bash
   git clone https://github.com/yourusername/HJ-DTA-DataBias.git
   cd HJ-DTA-DataBias
   ```

## Journal & Contact Info
[BASE: a web service for providing compound-protein binding affinity prediction datasets with reduced similarity bias](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-024-05968-3).

Hyojin Son*, Sechan Lee, Jaeuk Kim, Haangik Park, Myeong-Ha Hwang and Gwan-Su Yi†
- Journal: BMC Bioinformatics
- First Author(*): hyojin0912@kaist.ac.kr
- Corresponding Author(†): gwansuyi@kaist.ac.kr

## Acknowledgement
This work was supported by the BK-21 program through National Research Foundation of Korea (NRF) under Ministro of Education.

## License
The code in this repository is licensed under the MIT License. See the [LICENSE](./LICENSE) file for more details.

The data in the `data` folder is licensed under the CC BY 4.0 International license. See the [LICENSE](./LICENSE) file for more details.
