def gene_comb(buscos_folder, spicies_table):
    import pandas as pd 
    import os 
    all_types = []  # Список датафреймов
    type_names = []  # Список названий видов
    path_list = []  # Список путей до файлов

    # Добавляем все пути до файлов в указанной
    # директории
    for filename in os.scandir(buscos_folder):
        if filename.is_file() and filename.path.endswith(".tsv"):
            path_list.append(filename.path)

    table_sp = pd.read_csv(spicies_table, sep=';')
    gene_num = len(table_sp)

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
    return final

def find_all_missing():
    #Создаем серию в которой строки с комбинацией МММ имеют значение True
    ind_table = unated_table.eq('Missing', axis = 'columns').all(axis = 'columns')
    #Создаем список в котором лежат id генов отсутствующих у всех видов
    ind_list = unated_table.loc[ind_table].index