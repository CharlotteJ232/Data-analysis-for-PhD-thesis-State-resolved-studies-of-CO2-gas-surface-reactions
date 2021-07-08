import sys
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import \
    NavigationToolbar2QT as NavigationToolbar
import matplotlib.pyplot as plt
from matplotlib.figure import figaspect
import numpy as np
from PIL import Image
from matplotlib.ticker import NullFormatter


class Window(QWidget):
    def __init__(self, parent=None):
        super(Window, self).__init__(parent)

        self.figure = plt.figure()
        plt.gray()
        self.figure.tight_layout()
        self.canvas = FigureCanvas(self.figure)

        x = np.linspace(0, 2 * np.pi, 100)
        self.ori_data = (np.cos(x), np.sin(x))
        self.v = 180
        self.changed = 0
        self.has_rotated = 0
        self.toolbar = NavigationToolbar(self.canvas, self)

        # Just some button connected to `plot` method
        self.plot_button = QPushButton('Plot Image')
        self.plot_button.clicked.connect(self.plot_origin)
        self.plot_button.setDisabled(True)
        self.plot_check_box = QCheckBox('Plot origin?')

        self.open_button = QPushButton('Open Image')
        self.open_button.setStatusTip('open an image')
        self.open_button.clicked.connect(self.show_openfile)

        self.reverse_button = QPushButton('Reverse Image')
        self.reverse_button.setStatusTip('reverse the color')
        self.reverse_button.clicked.connect(self.reverse)
        self.reverse_button.setDisabled(True)

        self.hinst_button = QPushButton('Hinst Plot')
        self.hinst_button.clicked.connect(self.hinst)
        self.hinst_button.setDisabled(True)
        self.hinst_text = QLineEdit()

        self.hinst_text.setText(str(self.v))
        self.hinst_text.setStatusTip('larger value gets more points')
        self.hinst_refresh = QPushButton('Refresh')
        self.hinst_refresh.setStatusTip('update the zoomed image')
        self.hinst_refresh.setDisabled(True)
        self.hinst_refresh.clicked.connect(self.hinst_re)

        self.text = QTextEdit()

        self.save_img_button = QPushButton('Save Image')
        self.save_img_button.clicked.connect(self.save_img)
        self.save_data_button = QPushButton('Save Data')
        self.save_data_button.clicked.connect(self.save_data)

        self.save_data_button.setDisabled(True)
        self.save_img_button.setDisabled(True)

        mainbox = QHBoxLayout()

        vbox = QVBoxLayout()
        vbox.addWidget(self.toolbar)
        vbox.addWidget(self.canvas)

        hinst_group = QGroupBox('Hinst Plot')
        hinst_group.setFixedSize(250, 150)
        hinst_box = QGridLayout()
        hinst_box.addWidget(self.hinst_button, 0, 0)
        self.hinst_button.setStatusTip('Plot the full range of the hinst image')
        hinst_box.addWidget(QLabel('threhold: '), 1, 0)
        hinst_box.addWidget(self.hinst_refresh, 0, 1)
        hinst_box.addWidget(self.hinst_text, 1, 1)
        hinst_box.addWidget(self.save_data_button, 2, 0)
        hinst_group.setLayout(hinst_box)

        edge_group = QGroupBox('Find Edge')
        edge_group.setFixedSize(250, 150)
        edge_box = QGridLayout()
        edge_box.addWidget(
            QLabel('Edge line (left, top, width, height)'), 0, 0)
        self.edge_text = QLineEdit()
        self.edge_text.setStatusTip('the edge of the LEED screen')
        self.edge_text.setText('500, 500, 1500, 1500')
        edge_box.addWidget(self.edge_text, 1, 0)
        self.edge_plot_button = QPushButton('Plot edge')
        self.edge_plot_button.clicked.connect(self.edge_plot)
        edge_box.addWidget(self.edge_plot_button, 2, 0)
        self.edge_plot_button.setDisabled(True)
        self.edge_cut_button = QPushButton('Cut Edge')
        self.edge_cut_button.setDisabled(True)
        self.edge_cut_button.clicked.connect(self.edge_cut)
        edge_box.addWidget(self.edge_cut_button, 3, 0)
        edge_group.setLayout(edge_box)

        correct_group = QGroupBox('Correct Image')
        correct_group.setFixedSize(250, 120)
        correct_box = QGridLayout()
        correct_box.addWidget(QLabel('compress ratio: '), 0, 0)
        self.correct_ratio = QLineEdit()
        self.correct_ratio.setText('2')
        self.correct_ratio.setStatusTip('larger value makes the image smaller')
        correct_box.addWidget(self.correct_ratio, 1, 0)
        self.correct_button = QPushButton('Correct Image')
        self.correct_button.setStatusTip('correct the distorted image')
        self.correct_button.clicked.connect(self.correct)
        self.correct_button.setDisabled(True)
        correct_box.addWidget(self.correct_button, 2, 0)
        self.correct_bar = QProgressBar()
        self.correct_bar.setMaximum(100)
        self.correct_timer = QBasicTimer()
        self.correct_step = 0
        correct_box.addWidget(self.correct_bar, 3, 0)
        correct_group.setLayout(correct_box)

        rotate_group = QGroupBox('Rotate')
        rotate_group.setFixedSize(250, 100)
        rotate_box = QGridLayout()
        rotate_box.addWidget(QLabel('anti-clockwise'), 0, 0)
        self.rotate_text = QLineEdit()
        self.rotate_text.setText('0.0')
        rotate_box.addWidget(self.rotate_text, 1, 0)
        self.rotate_button = QPushButton('Rotate Image')
        self.rotate_button.clicked.connect(self.rotate)
        self.rotate_button.setDisabled(True)
        rotate_box.addWidget(self.rotate_button, 2, 0)
        rotate_group.setLayout(rotate_box)

        vbox2 = QVBoxLayout()
        vbox2.addWidget(self.open_button)
        vbox2.addWidget(self.plot_button)
        vbox2.addWidget(self.plot_check_box)
        vbox2.addWidget(self.reverse_button)
        vbox2.addWidget(edge_group)
        vbox2.addWidget(correct_group)
        vbox2.addWidget(rotate_group)
        vbox2.addWidget(hinst_group)
        vbox2.addWidget(QLabel('save DPI: '))
        self.dpi_text = QLineEdit()
        self.dpi_text.setStatusTip('larger value means higher quality')
        self.dpi_text.setText('200')
        self.dpi_text.setFixedSize(250, 20)
        vbox2.addWidget(self.dpi_text)
        vbox2.addWidget(self.save_img_button)

        mainbox.addLayout(vbox)
        mainbox.addLayout(vbox2)
        #
        self.setLayout(mainbox)
        self.setGeometry(200, 200, 1000, 800)

        # self.show()

    def edge_plot(self):

        self.plot()
        edge = self.edge_text.text().split(',')
        self.left, self.top, self.width, self.height = int(edge[0]), int(
            edge[1]), int(edge[2]), int(edge[3])
        plt.plot([self.left, self.left], [0, self.ori_lim[0]], 'r--')
        plt.plot([0, self.ori_lim[1]], [self.top, self.top], 'g--')
        plt.plot([self.left + self.width, self.left + self.width],
                 [0, self.ori_lim[0]], 'r--')
        plt.plot([0, self.ori_lim[1]],
                 [self.top + self.height, self.top + self.height], 'g--')

        self.canvas.draw()
        self.edge_cut_button.setEnabled(True)

    def edge_cut(self):
        self.figure.clear()
        data, lim = self.get_data()
        edge = self.edge_text.text().split(',')
        self.left, self.top, self.width, self.height = int(edge[0]), int(
            edge[1]), int(edge[2]), int(edge[3])
        ax = self.figure.add_subplot(111)
        ax.imshow(data[self.top:self.height + self.top,
                  self.left:self.left + self.width])
        self.canvas.draw()
        self.hinst_refresh.setDisabled(True)
        self.save_data_button.setDisabled(True)

    def plot(self):
        self.figure.clear()
        data, lim = self.get_data()
        ax = self.figure.add_subplot(111)
        ax.imshow(data)

        self.canvas.draw()
        self.hinst_refresh.setDisabled(True)
        self.save_data_button.setDisabled(True)

    def plot_origin(self):
        if self.plot_check_box.isChecked():
            self.figure.clear()
            data = self.ori_data
            lim = self.ori_lim
            ax = self.figure.add_subplot(111)
            ax.imshow(data)

            self.canvas.draw()
            self.hinst_refresh.setDisabled(True)
            self.save_data_button.setDisabled(True)
        else:
            self.plot()

    def reverse(self):
        if self.changed:
            self.new_data = 255 - self.new_data
        if self.has_rotated:
            self.rot_data = 255 - self.rot_data
        self.ori_data = 255 - self.ori_data
        self.plot()

    def get_data(self):
        pass

    def show_openfile(self):
        fn = QFileDialog.getOpenFileName()
        if fn[0]:
            self.ori_data = np.array(Image.open(fn[0]).convert('L'))
            self.ori_lim = self.ori_data.shape
            self.changed = 0
            self.has_rotated = 0
            self.plot_origin()
            self.setWindowTitle(fn[0])
            self.reverse_button.setEnabled(True)
            self.plot_button.setEnabled(True)
            self.hinst_button.setEnabled(True)
            self.save_img_button.setEnabled(True)
            self.correct_button.setEnabled(True)
            self.edge_plot_button.setEnabled(True)
            self.correct_button.setEnabled(True)
            self.rotate_button.setEnabled(True)

    def save_data(self):
        fn = QFileDialog.getSaveFileName()
        if fn[0]:
            l = max(len(self.h1[0]), len(self.h2[0]))
            with open(fn[0], 'w') as f:
                f.write('x_pos\tx_hinst\ty_pos\ty_hinst\n')
                for i in range(l):
                    if i < len(self.h1[0]):
                        f.write('%d\t' % self.h1[1][i])
                        f.write('%d\t' % self.h1[0][i])
                    else:
                        f.write('0\t0\t')
                    if i < len(self.h2[0]):
                        f.write('%d\t' % self.h2[1][i])
                        f.write('%d\n' % self.h2[0][i])
                    else:
                        f.write('0\t0\n')

    def hinst(self):
        self.v = int(self.hinst_text.text())
        v = self.v
        data, lim = self.get_data()

        try:
            sc1 = np.where(data < v)
        except Exception as e:
            print(e)
        self.x, self.y = sc1[1], -sc1[0]

        nullfmt = NullFormatter()  # no labels

        self.figure.clf()

        self.ax2 = self.figure.add_subplot(223)
        self.ax2.set_position([0.1, 0.1, 0.6, 0.6])
        self.ax2.scatter(self.x, self.y, s=1)
        self.ax2.set_xlim((0, lim[1]))
        self.ax2.set_ylim((-lim[0], 0))

        self.ax1 = self.figure.add_subplot(221)
        self.ax3 = self.figure.add_subplot(224)
        self.ax1.set_position([0.1, 0.8, 0.6, 0.15])
        self.ax1.set_yticklabels([])
        self.ax3.set_position([0.8, 0.1, 0.15, 0.6])
        self.ax3.set_xticklabels([])

        binwidth = 1

        bins1 = np.arange(0, lim[1] + binwidth, binwidth)
        bins2 = np.arange(-lim[0], 0 + binwidth, binwidth)
        self.h1, self.h2 = np.histogram(self.x, bins=bins1), \
                           np.histogram(self.y, bins=bins2)
        self.ax1.plot(self.h1[1][:-1], self.h1[0])
        self.ax3.plot(self.h2[0], self.h2[1][:-1])
        self.ax1.set_xlim(self.ax2.get_xlim())
        self.ax3.set_ylim(self.ax2.get_ylim())

        self.canvas.draw()
        self.hinst_refresh.setEnabled(True)
        self.save_data_button.setEnabled(True)

    def hinst_re(self):
        # print(self.ax2.get_ylim())
        self.v = int(self.hinst_text.text())
        data, lim = self.get_data()

        limx = self.ax2.get_xlim()
        limy = self.ax2.get_ylim()
        limx = np.array([int(limx[0]), int(limx[1])])
        limy = np.array([int(limy[0]), int(limy[1])])

        sv = np.where(data[-limy[1]:-limy[0], limx[0]:limx[1]] < self.v)
        self.x, self.y = sv[1] + limx[0], -sv[0] + limy[1]

        self.ax2.clear()
        self.ax1.clear()
        self.ax1.set_yticklabels([])
        self.ax3.clear()
        self.ax3.set_xticklabels([])

        self.ax2.scatter(self.x, self.y, s=1)
        self.ax2.set_xlim(limx)
        self.ax2.set_ylim(limy)

        binwidth = 1

        bins1 = np.arange(limx[0], limx[1] + binwidth, binwidth)
        bins2 = np.arange(limy[0], limy[1] + binwidth, binwidth)
        self.h1, self.h2 = np.histogram(self.x,
                                        bins=bins1), \
                           np.histogram(self.y, bins=bins2)

        self.ax1.plot(self.h1[1][:-1], self.h1[0])
        self.ax3.plot(self.h2[0], self.h2[1][:-1])

        self.ax1.set_xlim(limx)
        self.ax3.set_ylim(limy)

        self.canvas.draw()

    def save_img(self):
        fn = QFileDialog.getSaveFileName()
        dpi = int(self.dpi_text.text())
        if fn[0]:
            print(fn[0])
        plt.savefig(fn[0], dpi=dpi)

    def correct(self):
        if self.has_rotated:
            data = self.rot_data
            lim = self.rot_lim
        else:
            data = self.ori_data
            lim = self.ori_lim
        ratio = int(self.correct_ratio.text())
        im = data[self.top:self.height + self.top,
             self.left:self.left + self.width]
        im3 = np.array(
            Image.fromarray(im).resize(
                (int(im.shape[0] / ratio),
                 int(im.shape[1] / ratio))))
        xc, yc = im3.shape[1] / 2, im3.shape[0] / 2
        r, d, R = 4.75, 0.25, 8.3
        real_im = np.ones(im3.shape, dtype=im3.dtype) * 254
        r, d = r / R * sum(im3.shape) / 2, d / R * sum(im3.shape) / 2

        total = im3.shape[0] * im3.shape[1]
        for i in range(im3.shape[1]):
            QApplication.processEvents()
            for j in range(im3.shape[0]):
                QApplication.processEvents()
                real_i = r * (i - xc) / np.sqrt((r + d) ** 2 +
                                                (i - xc) ** 2 + (
                                                        j - yc) ** 2) + xc
                real_j = r * (j - yc) / np.sqrt((r + d) ** 2 +
                                                (i - xc) ** 2 + (
                                                        j - yc) ** 2) + yc

                real_im[int(real_j), int(real_i)] = im3[j, i]
                self.correct_step += 1
                if (self.correct_step * 100) % total == 0:
                    self.correct_bar.setValue(int(self.correct_step * 100 / total))

        self.new_data = real_im
        self.new_lim = self.new_data.shape
        self.changed = 1
        self.has_rotated = 0
        self.correct_step = 0
        self.plot()

    def rotate(self):
        if self.changed:
            data = self.new_data
            lim = self.new_lim
        else:
            data = self.ori_data
            lim = self.ori_lim
        angle = float(self.rotate_text.text())
        im = Image.fromarray(data)
        im2 = im.convert('RGBA')
        rot = im2.rotate(angle, expand=1)
        fff = Image.new('RGBA', rot.size, (255,) * 4)
        out = Image.composite(rot, fff, rot)

        self.rot_data = np.array(out.convert('L'))
        self.rot_lim = self.rot_data.shape
        self.has_rotated = 1
        self.plot()

    def get_data(self):
        if self.changed:
            if self.has_rotated:
                data = self.rot_data
                lim = self.rot_lim
            else:
                data = self.new_data
                lim = self.new_lim
        else:
            if self.has_rotated:
                data = self.rot_data
                lim = self.rot_lim
            else:
                data = self.ori_data
                lim = self.ori_lim
        return data, lim


class MainWindow(QMainWindow):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.plot_widget = Window(self)
        self.setCentralWidget(self.plot_widget)
        self.statusBar()

        self.setGeometry(200, 200, 1000, 800)

        openFile = QAction(QIcon('open.png'), 'Open', self)
        openFile.setShortcut('Ctrl+O')
        openFile.setStatusTip('open new File')
        openFile.triggered.connect(self.showDialog)

        self.plot_widget.plot_check_box.setStatusTip('Plot origin image?')

        # menubar = self.menuBar()
        # file_menu = menubar.addMenu('File')
        # file_menu.addAction(openFile)

    def showDialog(self):
        fname = QFileDialog.geself.topenFileName(self, 'Open file', '.')
        if fname[0]:
            self.plot_widget.text.setText(fname[0])
            im = np.array(Image.open(fname[0]).convert('L'))
            self.plot_widget.data = im
            self.plot_widget.plot()


if __name__ == '__main__':
    app = QApplication(sys.argv)

    try:
        main = MainWindow()
        main.show()
    except Exception as e:
        print(e)

    sys.exit(app.exec_())

    
    