# 📸 Image Restoration Project: Fixed Point Iteration (Richardson-Lucy)

This project implements a digital image deblurring system using the **Richardson-Lucy (RL)** algorithm. The method is treated as a **Fixed Point Problem**, directly applying the theory from **Chapter 3** of the Numerical Analysis course.

---

## 1. The Problem
In image processing, blur is modeled as a convolution:
$$y = Hx + n$$
Where **y** is the blurred image, **H** is the blur kernel (operator), and **x** is the original sharp image. Since we cannot simply "divide" by the blur, we use an iterative approach to recover the details.

## 2. Theoretical Solution (Chapter 3: Fixed Point)
We solve the equation by transforming it into a fixed point form:
$$x = g(x)$$

Using the Richardson-Lucy update rule:
$$x_{k+1} = x_k \cdot \left[ H^* \left( \frac{y}{Hx_k} \right) \right]$$

### ✅ Verification of Conditions (Theorem 4.1)
Before processing, the script automatically verifies the conditions required for convergence:
1. **Stability ($g([a, b]) \subseteq [a, b]$):** We verify the image and kernel are normalized. The script uses clamping to keep pixel values within the $[0, 1]$ interval.
2. **Contractivity:** The multiplicative update ensures the error decreases at each step, moving the estimate toward the unique fixed point.
3. **Normalization:** The sum of the kernel coefficients is checked (must equal 1) to ensure energy conservation.

---

## 3. How to Use
1. **Prepare your file:** Place your image in the root folder of the project.
2. **Default Image:** By default, the code looks for a file named `input.jpg`.
3. **Using PNG:** If you want to use a PNG:
   - Open `main.R`.
   - Change `input_file <- "input.jpg"` to `input_file <- "your_image.png"`.
   - The code will automatically handle the format.
4. **Run:** Execute the script in R. It will print the Condition Verification first, then process the image.

---

## 4. Project Outputs
The script generates the following files:
* **`blurred_light.jpg`**: The initial blurred version for reference.
* **`deblurred_final.jpg`**: The restored image after the Fixed Point iterations.
* **`comparison_final.jpg`**: A side-by-side view (Original vs. Blurred vs. Restored).
* **`rapport_final.txt`**: A full report including PSNR/SSIM scores and the mathematical verification logs.

---
**Course:** Numerical Analysis - ESI Sidi Bel Abbes
**Concept:** Resolution of Non-Linear Equations (Fixed Point)