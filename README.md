## Сборка
Для сборки проекта необходимо запустить ./build.sh предварительно прописав в него путь до vcpkg  
  
Данный скрипт соберет проект в директорию build, там будут находиться исполняемые файлы тестов tests и исполняемый файл для проведения симуляции работы метода порядковых статистик task-exec, он ожидает 7 аргументов: файл с порождающей матрицей, выходной файл, куда будут записаны результаты симуляции, количество итераций симуляции, шаг отношения сигнал/шум (в дб), его максимальное значение, 2 аргумента выбора алгоритма: 1 - использовать MakeCbtG, 0 - MakeCbtI; 1 - использовать CombCbtU, 0 - CombCbtV.

## Запуск
Для запуска тестов необходимо запустить run-tests.sh  
Для запуска симуляции и построения графиков run-task.sh  
