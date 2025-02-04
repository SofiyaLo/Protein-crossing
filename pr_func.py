import pandas as pd
import os
from matplotlib.colors import ListedColormap
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

def gene_comb(buscos_folder, spicies_table):
    all_types = []  # Список датафреймов
    type_names = []  # Список названий видов
    path_list = []  # Список путей до файлов

    # Добавляем все пути до файлов в указанной директории
    for filename in os.scandir(buscos_folder):
        if filename.is_file() and filename.path.endswith(".tsv"):
            path_list.append(filename.path)


    table_sp = pd.read_csv(spicies_table, sep=';') # Распаковываем таблицу из CSV формата
    gene_num = len(table_sp) # Создаем переменную с порядковым номером гена

    # Создаем список названий видов и список датафреймов полученных из busco
    for i in range(gene_num):
        type_name = table_sp.iloc[i, -1]
        file_name = table_sp.loc[table_sp.Название_вида == type_name, 'Название_файла'].iat[0]
        file = open(os.path.dirname(path_list[i] + '/' + file_name))
        print()
        df = pd.read_csv(file, sep='\t', index_col=[0], header=[2]). \
            replace({'Duplicated': 'Complete',
                     'Fragmented': 'Complete'})  # в качестве индекса датафреймов выступает название гена\
        # Заменяем для удобства подсчета
        type_names.append(type_name)
        all_types.append(df)

    # Удаляем дупликации
    for files in range(gene_num):
        all_types[files] = all_types[files][~all_types[files].index.duplicated(keep='first')]

    col_names = []  # Список названий столбцов, вид_статус

    for m in range(gene_num):
        column_name = (type_names[m] + '_status').capitalize()
        all_types[m] = all_types[m].rename(columns={'Status': column_name}).iloc[:, 0:1]
        # Отбираем только колонку со статусом гена
        col_names.append(column_name)

    unated_table = pd.concat(all_types, axis=1)  # Объедияем в одну таблицу колонки статусов

    final = unated_table.groupby(col_names).size().reset_index().rename(columns={0: 'freq'})
    print(final)
    return unated_table

def list_of_missing_genes(unated_table):
    # Создаем серию в которой строки с комбинацией МММ имеют значение True
    ind_table = unated_table.eq('Missing', axis='columns').all(axis='columns')

    # Создаем список в котором лежат id генов отсутствующих у всех видов
    ind_list = unated_table.loc[ind_table].index
    return ind_list

cmap = ListedColormap(['red', 'blue'])
def plot_status_heatmap(dataframe):
    # Map "Missing" на 0 и "Complete" на 1
    table_mapped = dataframe.map(lambda x: 0 if x == 'Missing' else 1 if x == 'Complete' else x)
    table_mapped = table_mapped.transpose()

    # Создаем тепловую карту
    #plt.figure(figsize=(10, 30))  # Adjust the figure size as needed
    sns.clustermap(table_mapped, col_cluster=False, cmap=cmap, cbar=True, method='ward', figsize=(8,30), cbar_pos=(0.91, 0.85, 0.03, 0.1))
    plt.setp(ax.yticklabels, rotation=45)

    # Меняем подписи цветовой шкалы
    colorbar = g.ax_heatmap.collections[0].colorbar
    colorbar.set_ticks([0, 1])  # Устанавливаем место расположения подписей относительно шкалы
    colorbar.set_ticklabels(['Missing', 'Complete'])  # Устанавливаем подписи

    # Подписываем оси графика
    ax.set_xlabel('Genes')
    ax.set_ylabel('Species')
    ax.set_title('Status Heatmap')
    plt.yticks(rotation=0)

    # Показать график
    plt.show()


def plot_clustered_heatmap(dataframe):

    # Map "Missing" на 0 и "Complete" на 1
    data_mapped = dataframe.map(lambda x: 0 if x == 'Missing' else 1 if x == 'Complete' else x)
    data_mapped = data_mapped.transpose()

    # Строим clustermap чтобы получить иерархическую кластеризацию для строк и столбцов
    g = sns.clustermap(data_mapped, cmap=cmap, cbar=True, cbar_pos=(0.8, 0.85, 0.03, 0.1), figsize=(13, 7))
    g.tick_params(axis='x', labelsize=6, rotation=45)

    # Меняем подписи цветовой шкалы
    colorbar = g.ax_heatmap.collections[0].colorbar
    colorbar.set_ticks([0, 1])  # Устанавливаем место расположения подписей относительно шкалы
    colorbar.set_ticklabels(['Missing', 'Complete'])  # Устанавливаем подписи

    # Показать clustermap
    plt.show()
    print(g.dendrogram_row.reordered_ind)


def PCA_graph(dataset, n_components=2):
    # Maппируем "Missing" на 0 и "Complete" на 1
    data_mapped = dataset.map(lambda x: 0 if x == 'Missing' else 1 if x == 'Complete' else x)
    data_mapped = data_mapped.transpose()
    # Введем индексы для видов
    data_mapped["spic_number"] = [i for i in range(1, len(data_mapped) + 1)]

    # Инициализируем PCA с нужным количеством компонент
    pca = PCA(n_components=n_components)

    # Обучаем PCA на маппированных нормализованных данных
    pca_result = pca.fit_transform(data_mapped.iloc[:, :-1])

    # Создаем датафрейм с результатами PCA
    pca_df = pd.DataFrame(pca_result, columns=[f'PC{i+1}' for i in range(n_components)])

    # Получим объясненную дисперсию для каждой из главных компонент
    explained_variance = pca.explained_variance_ratio_

    df = pd.DataFrame({"PC1": pca_result[:, 0],
                       "PC2": pca_result[:, 1],
                       "spic_number": data_mapped["spic_number"]})

    sns.scatterplot(x="PC1", y="PC2", hue="spic_number", data=df, alpha=0.75, palette="bright")
    #plt.show()

    return explained_variance