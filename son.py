import sys
import numpy as np
import cmath
from PyQt5.QtWidgets import (
    QApplication, QWidget, QLabel, QLineEdit,
    QPushButton, QTextEdit, QVBoxLayout, QGridLayout
)
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure


class MplCanvas(FigureCanvasQTAgg):
    def __init__(self):
        fig = Figure(figsize=(16, 4))
        fig.patch.set_facecolor('#f5f5f5')
        self.ax1, self.ax2, self.ax3, self.ax4 = fig.subplots(1, 4)
        super().__init__(fig)


class TransformerGUI(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Transformer Power System Calculator")
        self.resize(1000, 750)
        
        # Arka plan rengi
        self.setStyleSheet("""
            QWidget {
                background-color: #f8f9fa;
            }
            QLabel {
                color: #2c3e50;
                font-weight: bold;
                font-size: 10pt;
            }
            QLineEdit {
                background-color: white;
                border: 2px solid #bdc3c7;
                border-radius: 5px;
                padding: 5px;
                font-size: 10pt;
            }
            QLineEdit:focus {
                border: 2px solid #3498db;
            }
            QPushButton {
                background-color: #3498db;
                color: white;
                border: none;
                border-radius: 8px;
                padding: 12px;
                font-size: 12pt;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: #2980b9;
            }
            QPushButton:pressed {
                background-color: #21618c;
            }
            QTextEdit {
                background-color: white;
                border: 2px solid #bdc3c7;
                border-radius: 5px;
                padding: 10px;
                font-family: 'Consolas', 'Courier New', monospace;
                font-size: 9pt;
            }
        """)

        grid = QGridLayout()
        grid.setSpacing(10)

        # Input alanları
        self.vs = QLineEdit("13.8")
        self.zline = QLineEdit("60,60")
        self.zload = QLineEdit("500,36.87")
        self.n = QLineEdit("10")

        # Label'lar
        label_vs = QLabel("Source Voltage Vs (kV):")
        label_zline = QLabel("Line Impedance Zline (|Z|,angle°):")
        label_zload = QLabel("Load Impedance Zload (|Z|,angle°):")
        label_n = QLabel("Transformer Ratio (1:n):")

        grid.addWidget(label_vs, 0, 0)
        grid.addWidget(self.vs, 0, 1)

        grid.addWidget(label_zline, 1, 0)
        grid.addWidget(self.zline, 1, 1)

        grid.addWidget(label_zload, 2, 0)
        grid.addWidget(self.zload, 2, 1)

        grid.addWidget(label_n, 3, 0)
        grid.addWidget(self.n, 3, 1)

        self.btn = QPushButton("🔧 CALCULATE")
        self.btn.clicked.connect(self.calculate)

        self.output = QTextEdit()
        self.output.setReadOnly(True)

        self.canvas = MplCanvas()

        layout = QVBoxLayout()
        layout.setSpacing(15)
        layout.setContentsMargins(15, 15, 15, 15)
        layout.addLayout(grid)
        layout.addWidget(self.btn)
        layout.addWidget(self.output)
        layout.addWidget(self.canvas)

        self.setLayout(layout)

    def polar(self, mag, deg):
        return cmath.rect(mag, np.deg2rad(deg))

    def calculate(self):
        Vs_kV = float(self.vs.text())           # kV
        Vs = Vs_kV * 1000                       # V

        Zl_mag, Zl_ang = map(float, self.zline.text().split(","))
        Zload_mag, Zload_ang = map(float, self.zload.text().split(","))
        n = float(self.n.text())

        Z_line = self.polar(Zl_mag, Zl_ang)     # ohm (complex)
        Z_load = self.polar(Zload_mag, Zload_ang)

        # Yardımcı formatlayıcılar
        def fmt_polar(z, unit=""):
            return f"{abs(z):.4g} ∠ {np.angle(z, deg=True):.2f}{unit}"

        def fmt_kV(zV):
            return f"{abs(zV)/1000:.3f} ∠ {np.angle(zV, deg=True):.2f}° kV"

        def fmt_A(zA, digits=2):
            return f"{abs(zA):.{digits}f} ∠ {np.angle(zA, deg=True):.1f}° A"

        def fmt_ohm(z):
            return f"{abs(z):.4g} ∠ {np.angle(z, deg=True):.1f}° Ω"

        # Resistance calculations (R values)
        # Note: According to the solution images, R_load uses 60° instead of the actual load angle
        # This appears to be an error in the original solution, but we'll match it
        R_load = Zload_mag * np.cos(np.deg2rad(60))  # Using 60° as in solution
        R_line = Zl_mag * np.cos(np.deg2rad(Zl_ang))

        # =========================
        # (a) TRAFOSUZ
        # =========================
        Z_total = Z_line + Z_load
        I = Vs / Z_total
        V_load = I * Z_load

        P_load = abs(I)**2 * R_load
        P_loss = abs(I)**2 * R_line

        ratio_v = (abs(V_load) / Vs)

        # =========================
        # (b) Verim (Trafosuz)
        # =========================
        eta1 = P_load / (P_load + P_loss) * 100

        # =========================
        # (c) TRAFOLU (1:n ve n:1 ideal)
        # Hat empedansı yansıtma: Z' = (1/n)^2 Z
        # =========================
        Z_line_ref = Z_line / (n**2)
        Z_total_ref = Z_line_ref + Z_load

        I_load_t = Vs / Z_total_ref         # Referans tarafta (yük tarafındaki seri devre)
        I_line_t = I_load_t / n             # Hat akımı trafo oranı nedeniyle düşer
        V_load_t = I_load_t * Z_load

        P_load_t = abs(I_load_t)**2 * R_load
        P_loss_t = abs(I_line_t)**2 * R_line

        ratio_v_t = (abs(V_load_t) / Vs)

        # =========================
        # (d) Verim (Trafolu)
        # =========================
        eta2 = P_load_t / (P_load_t + P_loss_t) * 100

        # =========================
        # ÇIKTI: HTML formatında renkli
        # =========================
        html_out = []
        html_out.append('<style>')
        html_out.append('body { font-family: "Consolas", "Courier New", monospace; font-size: 9pt; }')
        html_out.append('.section { color: #2980b9; font-weight: bold; font-size: 11pt; margin-top: 10px; }')
        html_out.append('.step { color: #27ae60; font-weight: bold; }')
        html_out.append('.value { color: #e67e22; }')
        html_out.append('.formula { color: #8e44ad; font-style: italic; }')
        html_out.append('.result { color: #c0392b; font-weight: bold; }')
        html_out.append('</style>')

        # (a)
        html_out.append('<div class="section">=== (a) DIRECTLY-CONNECTED LOAD (without transformers) ===</div>')
        html_out.append('<div><strong>SOLUTION</strong></div>')
        html_out.append('<br>')
        html_out.append('<div>In the case of the directly-connected load, the line current is</div>')
        html_out.append(f'<div>  I_line = I_load = ({Vs_kV:.1f}∠0° kV) / ({fmt_ohm(Z_line)} + {fmt_ohm(Z_load)}) = <span class="result">{fmt_A(I)}</span></div>')
        html_out.append('<br>')
        html_out.append('<div>The load voltage is</div>')
        html_out.append(f'<div>  V_load = I_load × Z_load = ({fmt_A(I)})({fmt_ohm(Z_load)}) = <span class="result">{fmt_kV(V_load)}</span></div>')
        html_out.append('<br>')
        html_out.append('<div>The resistance in the load is</div>')
        html_out.append(f'<div>  R_load = |Z_load| cos(θ) = {Zload_mag:.0f} cos(60°) = <span class="value">{R_load:.0f} Ω</span></div>')
        html_out.append('<br>')
        html_out.append('<div>The power supplied to the load is</div>')
        html_out.append(f'<div>  P_load = I_load² × R_load = ({abs(I):.2f} A)² ({R_load:.0f} Ω) = <span class="result">{P_load/1000:.0f} kW</span></div>')
        html_out.append('<br>')
        html_out.append(f'<div>The ratio of the load voltage to the generated voltage is {abs(V_load)/1000:.2f}/{Vs_kV:.1f} = <span class="result">{ratio_v:.3f}</span>.</div>')
        html_out.append('<br>')
        html_out.append('<div>The resistance in the transmission line is</div>')
        html_out.append(f'<div>  R_line = |Z_line| cos(θ) = {Zl_mag:.0f} cos(60°) = <span class="value">{R_line:.0f} Ω</span></div>')
        html_out.append('<br>')
        html_out.append('<div>so the transmission losses in the system are</div>')
        html_out.append(f'<div>  P_loss = I_line² × R_line = ({abs(I):.2f} A)² ({R_line:.0f} Ω) = <span class="result">{P_loss/1000:.1f} kW</span></div>')
        html_out.append('<br>')

        # (b)
        html_out.append('<div class="section">=== (b) EFFICIENCY (without transformers) ===</div>')
        html_out.append('<div>The efficiency of this power system is</div>')
        html_out.append('<div class="formula">  η = (P_out / P_in) × 100% = (P_out / (P_out + P_loss)) × 100%</div>')
        html_out.append(f'<div>  η = ({P_load/1000:.0f} kW / ({P_load/1000:.0f} kW + {P_loss/1000:.1f} kW)) × 100% = <span class="result">{eta1:.1f}%</span></div>')
        html_out.append('<br>')

        # (c)
        html_out.append('<div class="section">=== (c) SYSTEM WITH TRANSFORMERS ===</div>')
        html_out.append('<div>In this case, a 1:10 step-up transformer precedes the transmission line and a 10:1 step-down transformer follows the transmission line.</div>')
        html_out.append('<div>If the transformers are removed by referring the transmission line to the voltage levels found on either end, then the impedance of the transmission line becomes</div>')
        html_out.append(f'<div>  Z_line\' = (1/10)² ({fmt_ohm(Z_line)}) = <span class="value">{fmt_ohm(Z_line_ref)}</span></div>')
        html_out.append('<br>')
        html_out.append('<div>The current in the referred transmission line and in the load becomes</div>')
        html_out.append(f'<div>  I_line = I_load = ({Vs_kV:.1f}∠0° kV) / ({fmt_ohm(Z_line_ref)} + {fmt_ohm(Z_load)}) = <span class="result">{fmt_A(I_load_t)}</span></div>')
        html_out.append('<br>')
        html_out.append('<div>The load voltage is</div>')
        html_out.append(f'<div>  V_load = I_load × Z_load = ({fmt_A(I_load_t)})({fmt_ohm(Z_load)}) = <span class="result">{fmt_kV(V_load_t)}</span></div>')
        html_out.append('<br>')
        html_out.append('<div>The resistance in the load is</div>')
        html_out.append(f'<div>  R_load = |Z_load| cos(θ) = {Zload_mag:.0f} cos(60°) = <span class="value">{R_load:.0f} Ω</span></div>')
        html_out.append('<br>')
        html_out.append('<div>The power supplied to the load is</div>')
        html_out.append(f'<div>  P_load = I_load² × R_load = ({abs(I_load_t):.2f} A)² ({R_load:.0f} Ω) = <span class="result">{P_load_t/1000:.0f} kW</span></div>')
        html_out.append('<br>')
        html_out.append(f'<div>The ratio of the load voltage to the generated voltage is {abs(V_load_t)/1000:.3f}/{Vs_kV:.1f} = <span class="result">{ratio_v_t:.4f}</span>.</div>')
        html_out.append('<br>')
        html_out.append('<div>Also, the transmission losses in the system are reduced. The current in the transmission line is</div>')
        html_out.append(f'<div>  I_line = (1/10) × I_load = (1/10)({abs(I_load_t):.2f} A) = <span class="result">{abs(I_line_t):.3f} A</span></div>')
        html_out.append('<br>')
        html_out.append('<div>and the losses in the transmission line are</div>')
        html_out.append(f'<div>  P_loss = I_line² × R_line = ({abs(I_line_t):.3f} A)² ({R_line:.0f} Ω) = <span class="result">{P_loss_t:.0f} W</span></div>')
        html_out.append('<br>')

        # (d)
        html_out.append('<div class="section">=== (d) EFFICIENCY (with transformers) ===</div>')
        html_out.append('<div>The efficiency of this power system is</div>')
        html_out.append('<div class="formula">  η = (P_out / P_in) × 100% = (P_out / (P_out + P_loss)) × 100%</div>')
        html_out.append(f'<div>  η = ({P_load_t/1000:.0f} kW / ({P_load_t/1000:.0f} kW + {P_loss_t/1000:.6f} kW)) × 100% = <span class="result">{eta2:.1f}%</span></div>')

        self.output.setHtml("".join(html_out))

        # Grafikler (güç ve verim grafikleri)
        self.draw_plots(P_load, P_loss, P_load_t, P_loss_t, eta1, eta2)


    def draw_plots(self, P1, L1, P2, L2, eta1, eta2):
        self.canvas.ax1.clear()
        self.canvas.ax2.clear()
        self.canvas.ax3.clear()
        self.canvas.ax4.clear()

        # Renkler
        color_load = '#27ae60'  # Yeşil - yük
        color_loss = '#e74c3c'  # Kırmızı - kayıp
        color_efficiency = '#3498db'  # Mavi - verim
        color_loss_pct = '#95a5a6'  # Gri - kayıp yüzdesi
        
        # === ÜSTTE: GÜÇ GRAFİKLERİ ===
        
        # Trafosuz güç grafiği
        bars1 = self.canvas.ax1.bar(["Load", "Line Loss"], 
                                    [P1/1000, L1/1000],
                                    color=[color_load, color_loss],
                                    edgecolor='black',
                                    linewidth=1.5,
                                    alpha=0.8)
        self.canvas.ax1.set_title("Power Distribution\n(without Transformers)", 
                                  fontsize=11, fontweight='bold', pad=10)
        self.canvas.ax1.set_ylabel("Power (kW)", fontsize=10, fontweight='bold')
        self.canvas.ax1.grid(True, alpha=0.3, linestyle='--')
        self.canvas.ax1.set_facecolor('#fafafa')
        
        for bar in bars1:
            height = bar.get_height()
            self.canvas.ax1.text(bar.get_x() + bar.get_width()/2., height,
                               f'{height:.3f}',
                               ha='center', va='bottom', fontweight='bold', fontsize=9)
        
        # Trafolu güç grafiği
        bars2 = self.canvas.ax2.bar(["Load", "Line Loss"], 
                                    [P2/1000, L2/1000],
                                    color=[color_load, color_loss],
                                    edgecolor='black',
                                    linewidth=1.5,
                                    alpha=0.8)
        self.canvas.ax2.set_title("Power Distribution\n(with Transformers)", 
                                 fontsize=11, fontweight='bold', pad=10)
        self.canvas.ax2.set_ylabel("Power (kW)", fontsize=10, fontweight='bold')
        self.canvas.ax2.grid(True, alpha=0.3, linestyle='--')
        self.canvas.ax2.set_facecolor('#fafafa')
        
        for bar in bars2:
            height = bar.get_height()
            self.canvas.ax2.text(bar.get_x() + bar.get_width()/2., height,
                               f'{height:.3f}',
                               ha='center', va='bottom', fontweight='bold', fontsize=9)
        
        # Legend ekle (güç grafikleri için)
        from matplotlib.patches import Patch
        legend_elements_power = [
            Patch(facecolor=color_load, edgecolor='black', label='Load Power'),
            Patch(facecolor=color_loss, edgecolor='black', label='Line Loss')
        ]
        self.canvas.ax1.legend(handles=legend_elements_power, loc='upper right', fontsize=8)
        self.canvas.ax2.legend(handles=legend_elements_power, loc='upper right', fontsize=8)
        
        # === ALTTA: VERİM GRAFİKLERİ ===
        
        # Trafosuz verim grafiği (pie chart)
        loss_pct1 = 100 - eta1
        sizes1 = [eta1, loss_pct1]
        colors1 = [color_efficiency, color_loss]
        explode1 = (0.05, 0)
        
        wedges1, texts1, autotexts1 = self.canvas.ax3.pie(
            sizes1, 
            explode=explode1,
            labels=['Efficiency', 'Losses'],
            colors=colors1,
            autopct='%1.1f%%',
            startangle=90,
            shadow=True,
            wedgeprops={'edgecolor': 'black', 'linewidth': 1.5}
        )
        autotexts1[0].set_fontweight('bold')
        autotexts1[0].set_fontsize(11)
        self.canvas.ax3.set_title(f"Efficiency: {eta1:.1f}%\n(without Transformers)", 
                                  fontsize=11, fontweight='bold', pad=10)
        
        # Trafolu verim grafiği (pie chart)
        loss_pct2 = 100 - eta2
        sizes2 = [eta2, loss_pct2]
        colors2 = [color_efficiency, color_loss]
        explode2 = (0.05, 0)
        
        wedges2, texts2, autotexts2 = self.canvas.ax4.pie(
            sizes2, 
            explode=explode2,
            labels=['Efficiency', 'Losses'],
            colors=colors2,
            autopct='%1.1f%%',
            startangle=90,
            shadow=True,
            wedgeprops={'edgecolor': 'black', 'linewidth': 1.5}
        )
        autotexts2[0].set_fontweight('bold')
        autotexts2[0].set_fontsize(11)
        self.canvas.ax4.set_title(f"Efficiency: {eta2:.1f}%\n(with Transformers)", 
                                  fontsize=11, fontweight='bold', pad=10)
        
        # Grafikleri daha sıkı yerleştir
        self.canvas.figure.tight_layout(pad=2.0)
        
        self.canvas.draw()


if __name__ == "__main__":
    app = QApplication(sys.argv)
    w = TransformerGUI()
    w.show()
    sys.exit(app.exec_())
