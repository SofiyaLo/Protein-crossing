import pandas as pd
import os
from matplotlib.colors import ListedColormap
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA


# Функция для чтения файла Protienortho
def ortho_reader(file):
    ortho_df = pd.read_csv(file, sep='\t')
    categorical_columns = ortho_df.columns[3:19]

    def replace_gene_values(value):
        if value == '*':
            return 0  # Заменить * на 0
        elif isinstance(value, str):
            return len(value.split(','))  # Количество генов (разделенных запятыми)
        return 0  # Если значение не строка, заменить на 0

    # Применение функции replace_gene_values ко всем значениям в столбцах с генами
    ortho_df[categorical_columns] = ortho_df[categorical_columns].applymap(replace_gene_values)
    return ortho_df

def PCA_graph(dataset, n_components=2):
    # Map "Missing" на 0 и "Complete" на 1
    data_mapped = dataset.map(lambda x: 0 if x == '*' else 1 if x == 'Complete' else x)
    data_mapped = data_mapped.transpose()
    data_mapped["spic_number"] = [i for i in range(1, len(data_mapped) + 1)]

    # Инициализируем PCA с нужным количеством компонент
    pca = PCA(n_components=n_components)

    # Обучаем PCA на маппированных нормализованных данных
    pca_result = pca.fit_transform(data_mapped.iloc[:, :-1])



