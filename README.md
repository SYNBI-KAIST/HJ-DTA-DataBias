# BASE: a web service for providing compound-protein binding affinity prediction datasets with reduced similarity bias

## Overview
Deep learning-based drug-target affinity (DTA) prediction models have shown high performance but suffer from dataset bias. Our study investigates this bias using comprehensive databases and demonstrates that compound-protein binding affinity can often be predicted using compound features alone, due to high similarity among target proteins. We developed bias-reduced datasets by decreasing protein similarity between training and test sets, which improved model performance and balanced feature importance.

We introduce the Binding Affinity Similarity Explorer (BASE) web service, which offers bias-reduced datasets and prediction results to aid in the development of generalized and robust DTA models. BASE is freely available at https://fundis.kaist.ac.kr/base.

![Figure9](https://github.com/user-attachments/assets/18336449-4f46-4ffd-9164-b80008e92d5a)

## Journal & Contact Info
[Korea Advanced Institute Science and Technology(KAIST)](https://kaist.ac.kr/en/)

[Department of Bio and Brain Engineering(BBE)](https://bioeng.kaist.ac.kr/)

[Synergistic Bioinformatics Laboratory](https://synbi.kaist.ac.kr/)

Hyojin Son*, Sechan Lee, Jaeuk Kim, Haangik Park, Myeong-Ha Hwang and Gwan-Su Yi†
- Journal: BMC Bioinformatics
- Status: Under Review
- First Author(*): hyojin0912@kaist.ac.kr
- Corresponding Author(†): gwansuyi@kaist.ac.kr

## Acknowledgement
This work was supported by the BK-21 program through National Research Foundation of Korea (NRF) under Ministro of Education.

## License
The code in this repository is licensed under the MIT License. See the [LICENSE](./LICENSE) file for more details.

The data in the `data` folder is licensed under the CC0 1.0 Universal (CC0 1.0) Public Domain Dedication. See the [LICENSE](./LICENSE) file for more details.
