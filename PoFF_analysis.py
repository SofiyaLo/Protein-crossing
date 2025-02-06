import pandas as pd
import os
from matplotlib.colors import ListedColormap
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA



def ortho_unweight_reader(file, threshold):
    """Функция читает файл результата ProteinOrtho и убирает неинформативные колонки.
        Каждое значение гена заменяет бинарно на наличие/отсутствие и фильтрует по заданному порогу
        threshold, обозначащем количество ортологов в семействе генов"""

    ortho_df = pd.read_csv(file, sep='\t').iloc[:, 3:19]
    ortho_df = ortho_df.map(lambda x: 0 if x == '*' else 1)
    ortho_df['zero_count'] = (ortho_df == 0).sum(axis=1)

    # Оставляем только строки, в которых количество нулей >= threshold
    filtered_df = ortho_df[ortho_df['zero_count'] <= threshold].drop(columns=['zero_count'])  # Удалим столбец с количеством нулей

    return filtered_df


def ortho_weight_reader(file, threshold):
    """Функция читает файл результата ProteinOrtho и убирает неинформативные колонки.
        Каждое значение гена заменяет количеством паралогов и фильтрует по заданному порогу
        threshold, обозначащем количество ортологов в семействе генов"""

    ortho_df = pd.read_csv(file, sep='\t').iloc[:, 3:19]

    ortho_df = ortho_df.map(lambda x: 0 if x == '*' else len(x.split(',')))
    ortho_df['zero_count'] = (ortho_df == 0).sum(axis=1)

    # Оставляем только строки, в которых количество нулей >= пороговому значению
    filtered_df = ortho_df[ortho_df['zero_count'] <= threshold].drop(columns=['zero_count'])  # Удалим столбец с количеством нулей

    return filtered_df

def PCA_graph_spic(dataset, n_components=2):
    """Функция принимает датафрейм из ProteinOrtho преобразованный функциями ortho_weight_reader
     или ortho_unweight_reader. Возвращает результат PCA по семействам генов включая график"""

    #Транспонируем датасет, чтобы гены стали строками а семейства колонками
    dataset = dataset.transpose()
    #Добавляем порядковый номер вида
    dataset["spic_number"] = [i for i in range(1, len(dataset) + 1)]

    # Инициализируем PCA с нужным количеством компонент
    pca = PCA(n_components=n_components)

    # Обучаем PCA на маппированных нормализованных данных
    pca_result = pca.fit_transform(dataset.iloc[:, :-1])

    df = pd.DataFrame({"PC1": pca_result[:, 0],
                       "PC2": pca_result[:, 1],
                       "spic_number": dataset["spic_number"]})

    sns.scatterplot(x="PC1", y="PC2", hue="spic_number", data=df, alpha=0.75, palette="bright")
    plt.show()

def PCA_graph_genes(dataset, n_components=2):
    """Функция принимает датафрейм из ProteinOrtho преобразованный функциями ortho_weight_reader
     или ortho_unweight_reader. Возвращает результат PCA по видам включая график"""

    #Добавляем порядковый номер семейства
    dataset["gene_number"] = [i for i in range(1, len(dataset) + 1)]

    # Инициализируем PCA с нужным количеством компонент
    pca = PCA(n_components=n_components)

    # Обучаем PCA на маппированных нормализованных данных
    pca_result = pca.fit_transform(dataset.iloc[:, :-1])

    df = pd.DataFrame({"PC1": pca_result[:, 0],
                       "PC2": pca_result[:, 1],
                       "gene_number": dataset["gene_number"]})

    sns.scatterplot(x="PC1", y="PC2", hue="gene_number", data=df, alpha=0.75, palette="bright")
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=6)

    plt.show()
