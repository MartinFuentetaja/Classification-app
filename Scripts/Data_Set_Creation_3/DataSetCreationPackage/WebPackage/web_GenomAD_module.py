#####################################################################################################################################################################################################################
#En este script vamos a tratar de descargar la info relativa a un gen de genomAD (V2 y V3_non_V2)

from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver import FirefoxOptions
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.common.exceptions import NoSuchElementException, TimeoutException

import os
import time
import sys

#####################################################################################################################################################################################################################

def download_genomAD(genomAD_url, Data_path, max_attempts=3, timeout=50):
    #Setting driver options
    opts = FirefoxOptions()
    opts.add_argument("--headless")
    opts.set_preference("browser.download.folderList", 2)
    opts.set_preference("browser.download.manager.showWhenStarting", False)
    opts.set_preference("browser.download.dir", Data_path)
    opts.set_preference("browser.helperApps.neverAsk.saveToDisk", "application/x-gzip")
    
    driver = webdriver.Firefox(options = opts)  # Reemplaza esto con el controlador adecuado para tu navegador
    driver.get(genomAD_url)
    driver.maximize_window()
    #Impose sleep so It has time to load all the page
    attempt = 1
    while attempt <= max_attempts:
        try:
            button = WebDriverWait(driver, timeout).until(
                EC.presence_of_element_located((By.XPATH, "/html/body/div[1]/div[3]/div[2]/div/div[7]/div[4]/div[2]/button[1]"))
            )
            button.click()
            break  # Salir del bucle si se hace clic con éxito
        except (TimeoutException, NoSuchElementException):
            print(f"Intento {attempt} fallido. No se encontró el elemento o tiempo de espera agotado.")
            attempt += 1
            if attempt <= max_attempts:
                print("Reintentando...")
            else:
                driver.quit()
                sys.exit("Attempts exhausted. Unable to download GenomAD file.")  # Salir del programa
    
    driver.quit()
    filename = max([os.path.join(Data_path, f) for f in os.listdir(Data_path)],key=os.path.getctime)
    return filename

def download_genomAD_V3_non_V2(genomAD_V3_url, Data_path):
    #Setting driver options
    opts = FirefoxOptions()
    opts.add_argument("--headless")
    opts.set_preference("browser.download.folderList", 2)
    opts.set_preference("browser.download.manager.showWhenStarting", False)
    opts.set_preference("browser.download.dir", Data_path)
    opts.set_preference("browser.helperApps.neverAsk.saveToDisk", "application/x-gzip")

    driver = webdriver.Firefox(options = opts)  # Reemplaza esto con el controlador adecuado para tu navegador
    driver.get(genomAD_V3_url)
    driver.maximize_window()
    #Impose sleep so It has time to load all the page
    time.sleep(15)
    button = driver.find_element(By.XPATH, "/html/body/div[1]/div[3]/div[2]/div/div[7]/div[4]/div[2]/button[1]")
    button.click()
    driver.quit()
    filename = max([os.path.join(Data_path, f) for f in os.listdir(Data_path)],key=os.path.getctime)
    return filename

def scraping_genomAD_info(genomAD_url):
    opts = FirefoxOptions()
    opts.add_argument("--headless")
    driver = webdriver.Firefox(options = opts)  # Reemplaza esto con el controlador adecuado para tu navegador
    driver.get(genomAD_url)
    time.sleep(1.5)
    info = driver.find_element(By.XPATH, "/html/body/div/div[3]").text
    driver.quit()
    return info
